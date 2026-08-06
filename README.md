# OmegaRG_E

A CLAS12 RG-E analysis program that reconstructs the **ω meson** via the decay chain

```
ω → π⁺ π⁻ π⁰
         └─→ γ γ
```

from `REC::Particle` banks in `.hipo` DST files, applying standard DIS event selection along the way. The program also produces a set of diagnostic histograms (multiplicities, per-particle kinematics, cutflow counters) so you can see exactly where events are being lost at each stage of the analysis.

---

## Requirements

- Java 17+ (uses `switch` expressions)
- [COATJAVA](https://github.com/JeffersonLab/clas12-offline-software) libraries on the classpath:
   - `org.jlab.clas.physics` (`LorentzVector`)
   - `org.jlab.jnp.hipo4` (`HipoReader`, `Bank`, `Event`, `SchemaFactory`)
   - `org.jlab.groot` (`H1F`, `TCanvas`, `F1D`, `DataFitter`)

## Running

```bash
javac *.java
java OmegaRG_E
```

The input directory is currently hardcoded in `OmegaRG_E.java`:

```java
String path = "C:\\Users\\mateo\\OneDrive\\Desktop\\summer\\HipoStuff";
```

Every `.hipo` file in that folder is processed. Update this path for your own setup.

---

## File overview

The program is split into four files, all in the same (default) package:

| File | Responsibility |
|---|---|
| `OmegaRG_E.java` | Main driver — file reading, particle selection, cuts, DIS variables, ω reconstruction |
| `Histograms.java` | Declares and constructs every `H1F` histogram used in the analysis |
| `KinematicsUtils.java` | Static helper methods shared across the analysis |
| `Plotting.java` | All `TCanvas` setup/drawing, kept separate from the event loop |

### `OmegaRG_E.java`

The entry point. High-level flow:

1. **File loop** — for each `.hipo` file in the input folder:
   - Opens it with `HipoReader` (see [File reading](#file-reading) below).
   - Grabs the `REC::Particle` schema/bank.
   - **Event loop** — for each event:
      1. Sorts particle rows into `electrons`, `gammas`, `pips`, `pims` lists by PID, keeping only rows with `status < 0` (COATJAVA's convention for "good" reconstructed particles — double check this against the CLAS12 note for your pass version).
      2. Fills raw multiplicity histograms (`hNElec`, `hNGamma`, `hNPip`, `hNPim`) **before** any topology cut, so you can see the raw distribution of particle counts per event.
      3. **Topology cut** — requires exactly 1 electron, 2 photons, 1 π⁺, 1 π⁻ (matches the `hipo-utils` filter string `"11,211,-211,22,22"`).
      4. Picks the **highest-momentum electron** as the scattered electron.
      5. **Vertex (target) cut** — keeps only electrons with `vz` inside a hardcoded window (`-13 < vz < -7`). See [Known limitations](#known-limitations) below.
      6. Fills p/θ/φ/vz histograms for every particle in the topology.
      7. Computes **DIS variables** (Q², W, y) from the electron's energy and scattering angle, and cuts on them (`Q² > 1.0`, `W > 2.0`, `y < 0.85`).
      8. Computes the **e⁻–γ opening angle** for both photons, fills `hEGammaAngle` (before the cut, so the full distribution — including the collinear peak — is visible), then cuts on `MIN_EGAMMA_ANGLE` to reject bremsstrahlung photons that are nearly collinear with the electron.
      9. Combines the two photons into a **π⁰ candidate**, fills `hPi0` with the diphoton mass, then cuts on the π⁰ mass window (`|m - 0.134| < 0.025` GeV).
      10. Combines π⁺, π⁻, and the π⁰ candidate into an **ω candidate**, fills `hOmega` with the invariant mass **only if it falls in `0.5 < m < 1.0` GeV**, and fills the ω's own p/θ/φ (from the summed 4-vector) plus a proxy vz (from the π⁺ track — see note in `Histograms.java`).
   - Closes the reader.
2. **Cutflow summary** — prints how many events survived each stage (total read → has electron → topology → vz → DIS → angle → π⁰ window → filled into `hOmega`). **This is the first thing to check if histograms come out empty** — it tells you exactly which cut is killing your statistics.
3. **π⁰ mass fit** — fits a Gaussian to `hPi0` between 0.10–0.17 GeV and prints the fitted mean/sigma.
4. **Plotting** — calls `Plotting.drawAll(h)` to open all canvases.

#### File reading

The current version reads files with `HipoReader`, one file at a time (open → drain all events → close, repeat). An earlier attempt used `HipoChain` to treat every `.hipo` file as one continuous stream, but `chain.add()` and `chain.close()` didn't resolve against the installed library version, so it was reverted to the confirmed-working `HipoReader` pattern. If you want to revisit `HipoChain`, check your IDE's autocomplete on the actual class to see what methods your jar version exposes before wiring it back in.

### `Histograms.java`

A single `Histograms` object bundles every `H1F` used by the analysis as public fields, constructed once at the start of `main()`. Grouped as:

- **Physics**: `hPi0`, `hOmega`
- **DIS kinematics**: `hQ2`, `hW`, `hy`
- **Multiplicities**: `hNGamma`, `hNPip`, `hNPim`, `hNElec`
- **Per-particle kinematics** (p, θ, φ, vz each): electron, π⁺, π⁻, γ
- **ω system kinematics**: `hOmegaP`, `hOmegaTh`, `hOmegaPh`, `hOmegaVz`
- **Angle diagnostic**: `hEGammaAngle`

Momentum-axis ranges scale with the beam energy passed into the constructor.

**Axis ranges were retuned against real output.** After an initial pass based on generic CLAS12 kinematics, the ranges were re-tuned a second time using actual canvas screenshots from a run of the code, so each histogram's axis now frames where the data actually falls rather than a generic guess:

| Histogram | Range         | Why                                                                                                                            |
|---|---------------|--------------------------------------------------------------------------------------------------------------------------------|
| `hPi0` | 0.05–0.20 GeV | subtle bump around 0.12–0.15 GeV on top of background                                                                          |
| `hOmega` | 0.7–.9 GeV    | matches the fill condition itself (`0.5 < m < 1.0`); the distribution rises all the way to the upper edge with.                |
| `hQ2` | 0.5–4.5 GeV²  | rises from ~1.0, peaks ~1.75–2.0, tapers by ~4.5                                                                               |
| `hW` | 1.0–4.5 GeV   | rises from ~1.0, peaks ~4.0, sharp kinematic cutoff near the beam-energy-limited W_max (~4.3)                                  |
| `hy` | 0.0–1.0       | rises smoothly across nearly the whole axis, peaking ~0.88–0.90                                                                |
| `hNGamma` | 0–8 (8 bins)  | small nonzero tail out to ~7–8                                                                                                 |
| `hNPip` | 0–5 (5 bins)  |                                                                                                                                |
| `hNPim` | 0–4 (4 bins)  |                                                                                                                                |
| `hNElec` | 0–3 (3 bins)  | essentially all events have exactly 1                                                                                          |
| `hElecP` | 1.5–8.5 GeV   | sharp onset ~1.6, tapers to ~0 by ~8.3                                                                                         |
| `hElecTh` | 5–32°         | peaks ~11–12°, tapers to ~0 by ~32°                                                                                            |
| `hElecVz` | -13 to -6     | sharp peak at ~-7.3                                                                                                            |
| `hPipP` / `hPimP` | 0–6 / 0–5 GeV | peak ~0.6–0.8 / ~0.4–0.6, tapering out                                                                                         |
| `hPipTh` / `hPimTh` | 5–70°         | long non-negligible tail out to 90° in the raw data                                                                            |
| `hPipVz` / `hPimVz` | -13 to -5     | dominant peak ~-7.5 to -7.7; a smaller secondary bump near -1.5 exists but isn't framed as the "main" peak                     |
| `hGamP` | 0–2.5 GeV     | sharp peak ~0.1–0.2, tapers to ~0 by ~2.5                                                                                      |
| `hGamTh` | 2–38°         | peaks ~12–13°, tapers to ~0 by ~38°                                                                                            |
| `hGamVz` | -13 to -6     | matches the charged-track vz peaks                                                                                             |
| `hOmegaP` | 0–7 GeV       | peaks ~2.5–3, tapers to ~0 by ~7                                                                                               |
| `hOmegaTh` | 0–40°         | peaks ~15–16°, tapers to ~0 by ~40°                                                                                            |
| `hOmegaVz` | -13 to -5     | matches π⁺ vz (used as the ω vz proxy)                                                                                         |
| `hEGammaAngle` | 0–0.7 rad     | shows both the near-zero collinear-background spike (removed by `MIN_EGAMMA_ANGLE`) and a broader genuine hump around ~0.5 rad |

All φ histograms (`hElecPh`, `hPipPh`, `hPimPh`, `hGamPh`, `hOmegaPh`) are left at the full -180 to 180° since they cover the whole azimuthal range and show real detector-sector structure worth seeing in full.

### `KinematicsUtils.java`

Static helpers used throughout the event loop:

| Method | Purpose |
|---|---|
| `getVector(Bank b, int row)` | Builds a `LorentzVector` from `px, py, pz` plus a PID-based mass assumption (γ: 0, π±: 0.13957 GeV, e⁻: 0.000511 GeV) |
| `angleBetween(LorentzVector a, LorentzVector b)` | Opening angle (radians) between two particles' 3-momenta |
| `fillKinematics(hP, hTh, hPh, hVz, v, vz)` | Fills the standard p/θ(deg)/φ(deg)/vz histogram set for one particle |
| `M_p()` | Proton mass constant (0.938272 GeV) |

### `Plotting.java`

All plotting logic lives here, separate from the event loop, so `OmegaRG_E.main()` just does analysis and then calls `Plotting.drawAll(h)` once at the end. Produces 8 `TCanvas` windows:

1. **DIS** — Q², W, y, π⁰ mass
2. **Omega** — ω invariant mass
3. **Multiplicities** — N(e⁻), N(γ), N(π⁺), N(π⁻)
4. **Electron_Kinematics** — e⁻ p/θ/φ/vz
5. **Pip_Kinematics** — π⁺ p/θ/φ/vz
6. **Pim_Kinematics** — π⁻ p/θ/φ/vz
7. **Gamma_Kinematics** — γ p/θ/φ/vz
8. **Omega_Kinematics_and_Angle** — ω p/θ/φ, e⁻–γ opening angle

---

## Known limitations / open TODOs

- **Visible peak is very subtle** The ω invariant mass distribution currently rises monotonically to the 1.0 GeV edge of the fill window with no turnover — this looks like combinatorial background dominating, not a resolved resonance. Widening the fill condition in `OmegaRG_E.java` (`if (m > 0.5 && m < 1.0)`) to see if it turns over past 1.0 GeV is a reasonable next step, along with revisiting the upstream cuts (π⁰ window, e⁻–γ angle, DIS cuts) to improve signal-to-background before concluding the ω isn't there.
- **π⁰ peak is subtle.** `hPi0` shows only a small bump over background around 0.12–0.15 GeV rather than a clean peak, which may make the Gaussian fit in `OmegaRG_E.java` unstable — worth checking the fit result (mean/sigma) against the raw histogram shape before trusting it.
- **e⁻–γ angle has two features, not clearly separated yet.** There's a spike near 0 rad (collinear background, removed by `MIN_EGAMMA_ANGLE = 0.02`) and a separate broad hump around 0.5 rad. It's not yet confirmed which of these — if either cleanly — corresponds to signal vs. background for this dataset.
- **Vertex cut is single-target only.** The `vz` window (`-13, -7`) currently selects one target. To handle both LD2 and a solid target correctly, look up each file's run number against the "Golden" tab of the RG-E run spreadsheet (column G) and choose the matching vz window per run, instead of hardcoding one window for everything.
- **`status < 0` PID-quality assumption.** This is COATJAVA's usual convention for good reconstructed particles, but should be double-checked against the CLAS12 note for whichever pass version produced your files.
- **ω vertex is a proxy.** The ω is a reconstructed composite (sum of π⁺, π⁻, π⁰ 4-vectors) with no directly measured vertex, so `hOmegaVz` uses the π⁺ track's vz as a stand-in under the assumption that the charged tracks share a common production vertex. A smaller secondary bump near vz ≈ -1.5 shows up in both `hPipVz` and `hPimVz` that isn't yet explained (second target contamination? a different vertex-finding path?) — worth a closer look.
- **Hardcoded input path.** `path` in `OmegaRG_E.java` points to a specific local folder and should be parameterized (command-line arg or config) for portability.
- **`HipoChain` unresolved.** File reading uses `HipoReader` per file rather than `HipoChain`, since the latter's API didn't match this project's installed library version (see [File reading](#file-reading)).

---

## Reading the cutflow output

If histograms come out empty, the cutflow block printed at the end of the run tells you where:

```
===== EVENT CUTFLOW =====
Total events read        = ...
Events with >=1 electron = ...
Pass exact topology cut  = ...
Pass vz (target) cut     = ...
Pass DIS cuts            = ...
Pass e-gamma angle cut   = ...
Pass pi0 mass window     = ...
Filled into hOmega       = ...
==========================
```

Each line should be less than or equal to the line above it. A sharp drop at any given line points to that specific cut being too tight (or a bug in that stage of the selection) for your dataset.