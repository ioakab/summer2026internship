import org.jlab.clas.physics.LorentzVector;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.SchemaFactory;

import org.jlab.groot.math.F1D;
import org.jlab.groot.fitter.DataFitter;

import java.io.File;
import java.util.ArrayList;

// NOTE: KinematicsUtils, Histograms, and Plotting are all in this same
// default (unnamed) package, so they're visible here with no import
// needed - just call them as KinematicsUtils.getVector(...), etc.

public class OmegaRG_E {

    // Beam energy (RG-E)
    static final double EBEAM = 10.6;

    // Minimum electron-photon opening angle (radians) to reject
    // collinear bremsstrahlung photons faking a pi0
    static final double MIN_EGAMMA_ANGLE = 0.02;

    public static void main(String[] args) {

        String path = "/Users/mateobisbal/Desktop/clas12-offline-software/HipoFold";

        File folder = new File(path);
        File[] files = folder.listFiles((d, name) -> name.endsWith(".hipo"));

        if (files == null || files.length == 0) {
            System.out.println("No files found.");
            return;
        }

        // All histograms now live in the Histograms class - see Histograms.java
        Histograms h = new Histograms(EBEAM);

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
        //
        // NOTE: an earlier version of this tried to use a HipoChain to
        // read all files as one continuous stream. That class's exact
        // method names (add/close/hasNext/etc.) vary between versions of
        // the hipo4 library and weren't resolving in this project's jar,
        // so this uses the confirmed-working HipoReader-per-file pattern
        // instead: open each file, drain its events, close it, move on.
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
                    //if (status >= 0) continue;

                    if (pid == 11)   electrons.add(i);
                    if (pid == 22)   gammas.add(i);
                    if (pid == 211)  pips.add(i);
                    if (pid == -211) pims.add(i);
                }
                //rec.show();

                // Multiplicity plots - filled on every event that has a
                // REC::Particle bank, before any topology cut
                h.hNGamma.fill(gammas.size());
                h.hNPip.fill(pips.size());
                h.hNPim.fill(pims.size());
                h.hNElec.fill(electrons.size());

                //System.out.println(gammas.size());
                // =========================
                // EXACT TOPOLOGY CUT: e-, pi+, pi-, 2 gamma
                // Matches the hipo-utils filter string "11,211,-211,22,22"
                // =========================
                if (electrons.size() >= 1) nHasElectron++;
                if (electrons.size() < 1) continue;


                // =========================
                // ELECTRON SELECTION (highest momentum)
                // =========================
                int bestE = electrons.get(0);
                double maxP = 0;

                for (int i : electrons) {
                    LorentzVector e = KinematicsUtils.getVector(rec, i);
                    double p = e.p();
                    if (p > maxP) {
                        maxP = p;
                        bestE = i;
                    }
                }

                LorentzVector electron = KinematicsUtils.getVector(rec, bestE);
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

                // =========================
                // DIS VARIABLES
                // =========================
                double E = EBEAM;
                double Eprime = electron.e();

                double Q2 = 4 * E * Eprime *
                        Math.pow(Math.sin(electron.theta() / 2.0), 2);

                double nu = E - Eprime;
                double W2 = KinematicsUtils.M_p() * KinematicsUtils.M_p() + 2 * KinematicsUtils.M_p() * nu - Q2;
                double W = Math.sqrt(Math.max(W2, 0));
                double y = nu / E;

                h.hQ2.fill(Q2);
                h.hW.fill(W);
                h.hy.fill(y);

                // DIS CUTS
                if (Q2 < 1.0) continue;
                if (W < 2.0) continue;
                if (y > 0.85) continue;

                nPassDIS++;


                if (gammas.size() != 2) continue;
                if (pips.size() != 1) continue;
                if (pims.size() != 1) continue;

                nPassTopology++;

                int ip = pips.get(0);
                int im = pims.get(0);
                int g1idx = gammas.get(0);
                int g2idx = gammas.get(1);

                LorentzVector pip = KinematicsUtils.getVector(rec, ip);
                LorentzVector pim = KinematicsUtils.getVector(rec, im);
                LorentzVector g1  = KinematicsUtils.getVector(rec, g1idx);
                LorentzVector g2  = KinematicsUtils.getVector(rec, g2idx);

                double vzPip = rec.getFloat("vz", ip);

                // Fill raw kinematics for every particle in this topology
                KinematicsUtils.fillKinematics(h.hElecP, h.hElecTh, h.hElecPh, h.hElecVz, electron, vzElectron);
                KinematicsUtils.fillKinematics(h.hPipP, h.hPipTh, h.hPipPh, h.hPipVz, pip, vzPip);
                KinematicsUtils.fillKinematics(h.hPimP, h.hPimTh, h.hPimPh, h.hPimVz, pim, rec.getFloat("vz", im));
                KinematicsUtils.fillKinematics(h.hGamP, h.hGamTh, h.hGamPh, h.hGamVz, g1, rec.getFloat("vz", g1idx));
                KinematicsUtils.fillKinematics(h.hGamP, h.hGamTh, h.hGamPh, h.hGamVz, g2, rec.getFloat("vz", g2idx));


                // =========================
                // ANGLE: e- vs each photon
                // Filled BEFORE the cut is applied so you can see the full
                // distribution (including the collinear peak you're cutting
                // away), then the cut rejects near-collinear pairs.
                // =========================
                double angle1 = KinematicsUtils.angleBetween(electron, g1);
                double angle2 = KinematicsUtils.angleBetween(electron, g2);

                h.hEGammaAngle.fill(angle1);
                h.hEGammaAngle.fill(angle2);

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
                h.hPi0.fill(mpi0);

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
                    h.hOmega.fill(m);

                    // omega system kinematics (p/theta/phi from the summed
                    // 4-vector; vz taken from the pi+ track as a proxy, see
                    // note in Histograms.java).
                    h.hOmegaP.fill(omega.p());
                    h.hOmegaTh.fill(Math.toDegrees(omega.theta()));
                    h.hOmegaPh.fill(Math.toDegrees(omega.phi()));
                    h.hOmegaVz.fill(vzPip);

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

        fit.setParameter(0, h.hPi0.getMax());
        fit.setParameter(1, 0.135);
        fit.setParameter(2, 0.01);

        DataFitter.fit(fit, h.hPi0, "Q");

        System.out.println("\n===== pi0 FIT =====");
        System.out.println("Mean  = " + fit.getParameter(1));
        System.out.println("Sigma = " + fit.getParameter(2));

        // =========================
        // PLOTS - see Plotting.java
        // =========================
        Plotting.drawAll(h);

        System.out.println("Done.");
    }
}