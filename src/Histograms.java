import org.jlab.groot.data.H1F;

/**
 * Holds every H1F histogram used by OmegaRG_E in one place, so the main
 * analysis class doesn't have ~30 lines of "new H1F(...)" at the top.
 * Construct one instance of this at the start of main() and pass it
 * around (or just reference its public fields directly).
 */
public class Histograms {

    static final double EBEAM = 10.6;

    // Beam energy is passed in so histogram ranges (e.g. momentum axes)
    // scale correctly if you ever run this at a different beam energy.
    public Histograms(double ebeam) {

        // Physics
        // hPi0 (DIS canvas, bottom-right): flat-ish with a subtle bump
        // around 0.12-0.15 GeV; narrowed to frame that bump.
        hPi0   = new H1F("hPi0", "M(#gamma#gamma)", 120, 0.05, 0.20);
        // hOmega (Omega canvas): the fill condition itself only keeps
        // 0.5 < m < 1.0, and the observed distribution rises all the way
        // to that upper edge with no turnover yet, so the full fill
        // window is already "the hump" - left as is. If you want to see
        // whether it turns over past 1.0, the fill condition in
        // OmegaRG_E.java (not just this axis) needs to widen too.
        hOmega = new H1F("hOmega", "M(#pi^{+}#pi^{-}#pi^{0})", 120, 0, 1.2);

        // DIS kinematics
        // hQ2 (DIS canvas, top-left): rises from ~1.0, peaks ~1.75-2.0,
        // tapers out by ~4.5.
        hQ2 = new H1F("hQ2", "Q^{2}", 100, 0.5, 4.5);
        // hW (DIS canvas, top-right): rises from ~1.0, peaks ~4.0, sharp
        // kinematic cutoff by ~4.3 (close to the W_max set by beam energy).
        hW  = new H1F("hW", "W", 100, 1.0, 4.5);
        // hy (DIS canvas, bottom-left): rises smoothly across nearly the
        // whole 0-1 axis, peaking ~0.88-0.90 - already fills the frame.
        hy  = new H1F("hy", "y", 100, 0.0, 1.0);

        // Multiplicities
        // (Multiplicities canvas): each distribution is dominated by one
        // or two low bins with almost nothing beyond - tightened per
        // particle to match what's actually populated.
        hNGamma = new H1F("hNGamma", "N(#gamma) per event", 8, 0, 8);
        hNPip   = new H1F("hNPip",   "N(#pi^{+}) per event", 5, 0, 5);
        hNPim   = new H1F("hNPim",   "N(#pi^{-}) per event", 4, 0, 4);
        hNElec  = new H1F("hNElec",  "N(e^{-}) per event", 3, 0, 3);

        // Electron kinematics (Electron_Kinematics canvas)
        // p: sharp onset ~1.6, broad plateau, tapers to ~0 by ~8.3.
        hElecP  = new H1F("hElecP",  "e^{-} p",     100, 1.5, 8.5);
        // theta: rises from ~6, peaks ~11-12, tapers to ~0 by ~32.
        hElecTh = new H1F("hElecTh", "e^{-} #theta",100, 5, 32);
        hElecPh = new H1F("hElecPh", "e^{-} #phi",  100, -180, 180);
        // vz: sharp peak at ~-7.3 with a small rising tail from -13.
        hElecVz = new H1F("hElecVz", "e^{-} vz",    100, -13, -6);

        // pi+ kinematics (Pip_Kinematics canvas)
        // p: peaks ~0.6-0.8, tapers to ~0 by ~6.
        hPipP  = new H1F("hPipP",  "#pi^{+} p",     100, 0, 6.0);
        // theta: rises from ~5, peaks ~15, long tail still nonzero at 90.
        hPipTh = new H1F("hPipTh", "#pi^{+} #theta",100, 5, 70);
        hPipPh = new H1F("hPipPh", "#pi^{+} #phi",  100, -180, 180);
        // vz: dominant peak at ~-7.7 (a smaller secondary bump near -1.5
        // is visible but this frames the main peak).
        hPipVz = new H1F("hPipVz", "#pi^{+} vz",    100, -13, -5);

        // pi- kinematics (Pim_Kinematics canvas)
        // p: peaks ~0.4-0.6, tapers to ~0 by ~5.
        hPimP  = new H1F("hPimP",  "#pi^{-} p",     100, 0, 5.0);
        // theta: rises from ~5, peaks ~23, long tail still nonzero at 90.
        hPimTh = new H1F("hPimTh", "#pi^{-} #theta",100, 5, 70);
        hPimPh = new H1F("hPimPh", "#pi^{-} #phi",  100, -180, 180);
        hPimVz = new H1F("hPimVz", "#pi^{-} vz",    100, -13, -5);

        // gamma kinematics (Gamma_Kinematics canvas, both photons fill these)
        // p: sharp peak ~0.1-0.2, tapers to ~0 by ~2.5.
        hGamP  = new H1F("hGamP",  "#gamma p",     100, 0, 2.5);
        // theta: rises from ~2, peaks ~12-13, tapers to ~0 by ~38.
        hGamTh = new H1F("hGamTh", "#gamma #theta",100, 2, 38);
        hGamPh = new H1F("hGamPh", "#gamma #phi",  100, -180, 180);
        // vz: sharp peak at ~-7.3, matches the charged-track vz peaks.
        hGamVz = new H1F("hGamVz", "#gamma vz",    100, -13, -6);

        // omega system kinematics (Omega_Kinematics_and_Angle canvas)
        // p: peaks ~2.5-3, tapers to ~0 by ~7.
        hOmegaP  = new H1F("hOmegaP",  "#omega p",     100, 0, 7.0);
        // theta: rises from 0, peaks ~15-16, tapers to ~0 by ~40.
        hOmegaTh = new H1F("hOmegaTh", "#omega #theta",100, 0, 40);
        hOmegaPh = new H1F("hOmegaPh", "#omega #phi",  100, -180, 180);
        hOmegaVz = new H1F("hOmegaVz", "#omega vz (from #pi^{+} track)", 100, -13, -5);

        // e- to photon opening angle (Omega_Kinematics_and_Angle canvas,
        // bottom-right). Two features: a spike near 0 (collinear
        // background the MIN_EGAMMA_ANGLE cut removes) and a broader
        // genuine hump centered ~0.45-0.55 that tapers out by ~0.7-0.8.
        // Framed to show both while trimming the empty tail past ~0.7.
        hEGammaAngle = new H1F("hEGammaAngle", "e^{-}-#gamma opening angle", 100, 0.0, 0.7);
    }

    public final H1F hPi0, hOmega;
    public final H1F hQ2, hW, hy;
    public final H1F hNGamma, hNPip, hNPim, hNElec;
    public final H1F hElecP, hElecTh, hElecPh, hElecVz;
    public final H1F hPipP, hPipTh, hPipPh, hPipVz;
    public final H1F hPimP, hPimTh, hPimPh, hPimVz;
    public final H1F hGamP, hGamTh, hGamPh, hGamVz;
    public final H1F hOmegaP, hOmegaTh, hOmegaPh, hOmegaVz;
    public final H1F hEGammaAngle;
}