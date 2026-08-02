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
        hPi0   = new H1F("hPi0", "M(#gamma#gamma)", 120, 0.0, 0.30);
        hOmega = new H1F("hOmega", "M(#pi^{+}#pi^{-}#pi^{0})", 120, 0.5, 1.0);

        // DIS kinematics
        hQ2 = new H1F("hQ2", "Q^{2}", 100, 0.0, 5.0);
        hW  = new H1F("hW", "W", 100, 1.0, 8.0);
        hy  = new H1F("hy", "y", 100, 0.0, 1.0);

        // Multiplicities
        hNGamma = new H1F("hNGamma", "N(#gamma) per event", 10, 0, 10);
        hNPip   = new H1F("hNPip",   "N(#pi^{+}) per event", 10, 0, 10);
        hNPim   = new H1F("hNPim",   "N(#pi^{-}) per event", 10, 0, 10);
        hNElec  = new H1F("hNElec",  "N(e^{-}) per event", 10, 0, 10);

        // Electron kinematics
        hElecP  = new H1F("hElecP",  "e^{-} p",     100, 0, ebeam);
        hElecTh = new H1F("hElecTh", "e^{-} #theta",100, 0, 50);
        hElecPh = new H1F("hElecPh", "e^{-} #phi",  100, -180, 180);
        hElecVz = new H1F("hElecVz", "e^{-} vz",    100, -20, 10);

        // pi+ kinematics
        hPipP  = new H1F("hPipP",  "#pi^{+} p",     100, 0, ebeam);
        hPipTh = new H1F("hPipTh", "#pi^{+} #theta",100, 0, 90);
        hPipPh = new H1F("hPipPh", "#pi^{+} #phi",  100, -180, 180);
        hPipVz = new H1F("hPipVz", "#pi^{+} vz",    100, -20, 10);

        // pi- kinematics
        hPimP  = new H1F("hPimP",  "#pi^{-} p",     100, 0, ebeam);
        hPimTh = new H1F("hPimTh", "#pi^{-} #theta",100, 0, 90);
        hPimPh = new H1F("hPimPh", "#pi^{-} #phi",  100, -180, 180);
        hPimVz = new H1F("hPimVz", "#pi^{-} vz",    100, -20, 10);

        // gamma kinematics (both photons fill the same histograms)
        hGamP  = new H1F("hGamP",  "#gamma p",     100, 0, ebeam);
        hGamTh = new H1F("hGamTh", "#gamma #theta",100, 0, 50);
        hGamPh = new H1F("hGamPh", "#gamma #phi",  100, -180, 180);
        hGamVz = new H1F("hGamVz", "#gamma vz",    100, -20, 10);

        // omega system kinematics (see note in OmegaRG_E about vz proxy)
        hOmegaP  = new H1F("hOmegaP",  "#omega p",     100, 0, ebeam);
        hOmegaTh = new H1F("hOmegaTh", "#omega #theta",100, 0, 50);
        hOmegaPh = new H1F("hOmegaPh", "#omega #phi",  100, -180, 180);
        hOmegaVz = new H1F("hOmegaVz", "#omega vz (from #pi^{+} track)", 100, -20, 10);

        // e- to photon opening angle (filled for both photons, before cut)
        hEGammaAngle = new H1F("hEGammaAngle", "e^{-}-#gamma opening angle", 100, 0.0, 1.0);
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