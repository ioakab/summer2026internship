import org.jlab.clas.physics.LorentzVector;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.SchemaFactory;

import org.jlab.groot.data.H1F;
import org.jlab.groot.ui.TCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.fitter.DataFitter;

import java.io.File;
import java.util.ArrayList;

public class OmegaRG_E {

    // Beam energy (RG-E)
    static final double EBEAM = 10.6;

    // Minimum electron-photon opening angle (radians) to reject
    // collinear bremsstrahlung photons faking a pi0
    static final double MIN_EGAMMA_ANGLE = 0.02;

    public static void main(String[] args) {

        String path = "C:\\Users\\mateo\\OneDrive\\Desktop\\summer\\HipoStuff";

        File folder = new File(path);
        File[] files = folder.listFiles((d, name) -> name.endsWith(".hipo"));

        if (files == null || files.length == 0) {
            System.out.println("No files found.");
            return;
        }

        // =========================
        // HISTOGRAMS - physics
        // =========================
        H1F hPi0   = new H1F("hPi0", "M(#gamma#gamma)", 120, 0.0, 0.30);
        H1F hOmega = new H1F("hOmega", "M(#pi^{+}#pi^{-}#pi^{0})", 120, 0.5, 1.0);

        // DIS kinematics
        H1F hQ2 = new H1F("hQ2", "Q^{2}", 100, 0.0, 5.0);
        H1F hW  = new H1F("hW", "W", 100, 1.0, 3.0);
        H1F hy  = new H1F("hy", "y", 100, 0.0, 1.0);

        // Multiplicities (filled BEFORE the exact-topology cut, so you can
        // see the raw distribution of how many of each particle show up)
        H1F hNGamma = new H1F("hNGamma", "N(#gamma) per event", 10, 0, 10);
        H1F hNPip   = new H1F("hNPip",   "N(#pi^{+}) per event", 10, 0, 10);
        H1F hNPim   = new H1F("hNPim",   "N(#pi^{-}) per event", 10, 0, 10);
        H1F hNElec  = new H1F("hNElec",  "N(e^{-}) per event", 10, 0, 10);

        // Per-particle kinematics: momentum, theta, phi, vz
        H1F hElecP  = new H1F("hElecP",  "e^{-} p",     100, 0, EBEAM);
        H1F hElecTh = new H1F("hElecTh", "e^{-} #theta",100, 0, 50);
        H1F hElecPh = new H1F("hElecPh", "e^{-} #phi",  100, -180, 180);
        H1F hElecVz = new H1F("hElecVz", "e^{-} vz",    100, -20, 10);

        H1F hPipP  = new H1F("hPipP",  "#pi^{+} p",     100, 0, EBEAM);
        H1F hPipTh = new H1F("hPipTh", "#pi^{+} #theta",100, 0, 90);
        H1F hPipPh = new H1F("hPipPh", "#pi^{+} #phi",  100, -180, 180);
        H1F hPipVz = new H1F("hPipVz", "#pi^{+} vz",    100, -20, 10);

        H1F hPimP  = new H1F("hPimP",  "#pi^{-} p",     100, 0, EBEAM);
        H1F hPimTh = new H1F("hPimTh", "#pi^{-} #theta",100, 0, 90);
        H1F hPimPh = new H1F("hPimPh", "#pi^{-} #phi",  100, -180, 180);
        H1F hPimVz = new H1F("hPimVz", "#pi^{-} vz",    100, -20, 10);

        H1F hGamP  = new H1F("hGamP",  "#gamma p",     100, 0, EBEAM);
        H1F hGamTh = new H1F("hGamTh", "#gamma #theta",100, 0, 50);
        H1F hGamPh = new H1F("hGamPh", "#gamma #phi",  100, -180, 180);
        H1F hGamVz = new H1F("hGamVz", "#gamma vz",    100, -20, 10);

        // =========================
        // DIAGNOSTIC COUNTERS
        // Use these to see exactly which cut is killing your events
        // =========================
        long nEventsTotal   = 0;
        long nHasElectron   = 0;
        long nPassTopology  = 0;
        long nPassVz        = 0;
        long nPassDIS       = 0;
        long nPassAngle     = 0;
        long nPassPi0       = 0;
        long nPassOmega     = 0;

        // =========================
        // FILE LOOP
        // =========================
        for (File file : files) {

            System.out.println("Processing: " + file.getName());

            HipoReader reader = new HipoReader();
            reader.open(file.getAbsolutePath());

            Event event = new Event();
            SchemaFactory factory = reader.getSchemaFactory();
            Bank rec = new Bank(factory.getSchema("REC::Particle"));

            // =========================
            // EVENT LOOP
            // =========================
            while (reader.hasNext()) {

                reader.nextEvent(event);
                event.read(rec);

                nEventsTotal++;

                int n = rec.getRows();
                if (n < 2) continue;

                ArrayList<Integer> gammas    = new ArrayList<>();
                ArrayList<Integer> pips      = new ArrayList<>();
                ArrayList<Integer> pims      = new ArrayList<>();
                ArrayList<Integer> electrons = new ArrayList<>();

                // =========================
                // PARTICLE SELECTION
                // =========================
                for (int i = 0; i < n; i++) {

                    int pid = rec.getInt("pid", i);
                    int status = rec.getInt("status", i);

                    // NOTE: this assumes status < 0 marks "good" reconstructed
                    // particles per the COATJAVA REC::Particle convention.
                    // Double check this against the CLAS12 note for your pass version.
                    if (status >= 0) continue;

                    if (pid == 11)   electrons.add(i);
                    if (pid == 22)   gammas.add(i);
                    if (pid == 211)  pips.add(i);
                    if (pid == -211) pims.add(i);
                }

                // Multiplicity plots - filled on every event that has a
                // REC::Particle bank, before any topology cut
                hNGamma.fill(gammas.size());
                hNPip.fill(pips.size());
                hNPim.fill(pims.size());
                hNElec.fill(electrons.size());

                // =========================
                // EXACT TOPOLOGY CUT: e-, pi+, pi-, 2 gamma
                // Matches the hipo-utils filter string "11,211,-211,22,22"
                // =========================
                if (electrons.size() >= 1) nHasElectron++;

                if (electrons.size() < 1) continue;
                if (gammas.size() != 2) continue;
                if (pips.size() != 1) continue;
                if (pims.size() != 1) continue;

                nPassTopology++;

                // =========================
                // ELECTRON SELECTION (highest momentum)
                // =========================
                int bestE = electrons.get(0);
                double maxP = 0;

                for (int i : electrons) {
                    LorentzVector e = getVector(rec, i);
                    double p = e.p();
                    if (p > maxP) {
                        maxP = p;
                        bestE = i;
                    }
                }

                LorentzVector electron = getVector(rec, bestE);
                double vzElectron = rec.getFloat("vz", bestE);

                // =========================
                // VERTEX CUT (TARGET SEPARATION)
                // TODO: this window currently selects ONE target only.
                // To handle both LD2 and the solid target, look up the
                // run number for this file against the "Golden" tab of the
                // RG-E run spreadsheet (column G) and pick the vz window
                // that matches that run's target instead of hardcoding it.
                // =========================
                if (vzElectron < -13 || vzElectron > -7) continue;

                nPassVz++;

                int ip = pips.get(0);
                int im = pims.get(0);
                int g1idx = gammas.get(0);
                int g2idx = gammas.get(1);

                LorentzVector pip = getVector(rec, ip);
                LorentzVector pim = getVector(rec, im);
                LorentzVector g1  = getVector(rec, g1idx);
                LorentzVector g2  = getVector(rec, g2idx);

                // Fill raw kinematics for every particle in this topology
                fillKinematics(hElecP, hElecTh, hElecPh, hElecVz, electron, vzElectron);
                fillKinematics(hPipP, hPipTh, hPipPh, hPipVz, pip, rec.getFloat("vz", ip));
                fillKinematics(hPimP, hPimTh, hPimPh, hPimVz, pim, rec.getFloat("vz", im));
                fillKinematics(hGamP, hGamTh, hGamPh, hGamVz, g1, rec.getFloat("vz", g1idx));
                fillKinematics(hGamP, hGamTh, hGamPh, hGamVz, g2, rec.getFloat("vz", g2idx));

                // =========================
                // DIS VARIABLES
                // =========================
                double E = EBEAM;
                double Eprime = electron.e();

                double Q2 = 4 * E * Eprime *
                        Math.pow(Math.sin(electron.theta() / 2.0), 2);

                double nu = E - Eprime;
                double W2 = M_p() * M_p() + 2 * M_p() * nu - Q2;
                double W = Math.sqrt(Math.max(W2, 0));
                double y = nu / E;

                hQ2.fill(Q2);
                hW.fill(W);
                hy.fill(y);

                // DIS CUTS
                if (Q2 < 1.0) continue;
                if (W < 2.0) continue;
                if (y > 0.85) continue;

                nPassDIS++;

                // =========================
                // ANGLE CUT: e- vs each photon
                // Reject events where either photon is nearly collinear
                // with the electron (radiative / Moller-like background)
                // =========================
                double angle1 = angleBetween(electron, g1);
                double angle2 = angleBetween(electron, g2);
                if (angle1 < MIN_EGAMMA_ANGLE) continue;
                if (angle2 < MIN_EGAMMA_ANGLE) continue;

                nPassAngle++;

                // =========================
                // gamma gamma -> pi0
                // =========================
                LorentzVector pi0 = new LorentzVector();
                pi0.add(g1);
                pi0.add(g2);

                double mpi0 = pi0.mass();
                hPi0.fill(mpi0);

                // pi0 mass window selection
                if (Math.abs(mpi0 - 0.134) > 0.025) continue;

                nPassPi0++;

                // =========================
                // omega -> pi+ pi- pi0   (filled ONCE per event)
                // =========================
                LorentzVector omega = new LorentzVector();
                omega.add(pip);
                omega.add(pim);
                omega.add(pi0);

                double m = omega.mass();
                if (m > 0.5 && m < 1.0) {
                    hOmega.fill(m);
                    nPassOmega++;
                }
            }

            reader.close();
        }

        // =========================
        // DIAGNOSTIC SUMMARY - read this first if histograms are blank
        // =========================
        System.out.println("\n===== EVENT CUTFLOW =====");
        System.out.println("Total events read        = " + nEventsTotal);
        System.out.println("Events with >=1 electron = " + nHasElectron);
        System.out.println("Pass exact topology cut  = " + nPassTopology);
        System.out.println("Pass vz (target) cut     = " + nPassVz);
        System.out.println("Pass DIS cuts            = " + nPassDIS);
        System.out.println("Pass e-gamma angle cut   = " + nPassAngle);
        System.out.println("Pass pi0 mass window     = " + nPassPi0);
        System.out.println("Filled into hOmega       = " + nPassOmega);
        System.out.println("==========================\n");

        // =========================
        // pi0 FIT
        // =========================
        F1D fit = new F1D("fit", "[amp]*gaus(x,[mean],[sigma])", 0.10, 0.17);

        fit.setParameter(0, hPi0.getMax());
        fit.setParameter(1, 0.135);
        fit.setParameter(2, 0.01);

        DataFitter.fit(fit, hPi0, "Q");

        System.out.println("\n===== pi0 FIT =====");
        System.out.println("Mean  = " + fit.getParameter(1));
        System.out.println("Sigma = " + fit.getParameter(2));

        // =========================
        // PLOTS
        // =========================
        TCanvas c1 = new TCanvas("DIS", 900, 600);
        c1.divide(2, 2);
        c1.cd(0); c1.draw(hQ2);
        c1.cd(1); c1.draw(hW);
        c1.cd(2); c1.draw(hy);
        c1.cd(3); c1.draw(hPi0);

        TCanvas c2 = new TCanvas("Omega", 800, 600);
        c2.draw(hOmega);

        TCanvas c3 = new TCanvas("Multiplicities", 900, 600);
        c3.divide(2, 2);
        c3.cd(0); c3.draw(hNElec);
        c3.cd(1); c3.draw(hNGamma);
        c3.cd(2); c3.draw(hNPip);
        c3.cd(3); c3.draw(hNPim);

        TCanvas c4 = new TCanvas("Electron_Kinematics", 900, 600);
        c4.divide(2, 2);
        c4.cd(0); c4.draw(hElecP);
        c4.cd(1); c4.draw(hElecTh);
        c4.cd(2); c4.draw(hElecPh);
        c4.cd(3); c4.draw(hElecVz);

        TCanvas c5 = new TCanvas("Pip_Kinematics", 900, 600);
        c5.divide(2, 2);
        c5.cd(0); c5.draw(hPipP);
        c5.cd(1); c5.draw(hPipTh);
        c5.cd(2); c5.draw(hPipPh);
        c5.cd(3); c5.draw(hPipVz);

        TCanvas c6 = new TCanvas("Pim_Kinematics", 900, 600);
        c6.divide(2, 2);
        c6.cd(0); c6.draw(hPimP);
        c6.cd(1); c6.draw(hPimTh);
        c6.cd(2); c6.draw(hPimPh);
        c6.cd(3); c6.draw(hPimVz);

        TCanvas c7 = new TCanvas("Gamma_Kinematics", 900, 600);
        c7.divide(2, 2);
        c7.cd(0); c7.draw(hGamP);
        c7.cd(1); c7.draw(hGamTh);
        c7.cd(2); c7.draw(hGamPh);
        c7.cd(3); c7.draw(hGamVz);

        System.out.println("Done.");
    }

    // =========================
    // Helper: fill p, theta, phi, vz histograms for one particle
    // =========================
    public static void fillKinematics(H1F hP, H1F hTh, H1F hPh, H1F hVz,
                                      LorentzVector v, double vz) {
        hP.fill(v.p());
        hTh.fill(Math.toDegrees(v.theta()));
        hPh.fill(Math.toDegrees(v.phi()));
        hVz.fill(vz);
    }

    // =========================
    // Helper: opening angle between two particles (radians)
    // =========================
    public static double angleBetween(LorentzVector a, LorentzVector b) {
        return Math.acos(
                a.vect().dot(b.vect()) / (a.vect().mag() * b.vect().mag())
        );
    }

    // =========================
    // LORENTZ VECTOR BUILDER
    // =========================
    public static LorentzVector getVector(Bank b, int row) {

        double px = b.getFloat("px", row);
        double py = b.getFloat("py", row);
        double pz = b.getFloat("pz", row);

        int pid = b.getInt("pid", row);

        double mass = switch (pid) {
            case 22 -> 0.0;
            case 211, -211 -> 0.13957;
            case 11 -> 0.000511;
            default -> 0.0;
        };

        LorentzVector v = new LorentzVector();
        v.setPxPyPzM(px, py, pz, mass);

        return v;
    }

    // proton mass
    public static double M_p() {
        return 0.938272;
    }
}