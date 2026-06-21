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

    public static void main(String[] args) {

        String path = "/home/teo/Documents/HipoFold";

        File folder = new File(path);
        File[] files = folder.listFiles((d, name) -> name.endsWith(".hipo"));

        if (files == null || files.length == 0) {
            System.out.println("No files found.");
            return;
        }

        // =========================
        // HISTOGRAMS
        // =========================
        H1F hPi0 = new H1F("hPi0", "M(#gamma#gamma)", 120, 0.0, 0.30);
        H1F hOmega = new H1F("hOmega", "M(#pi^{+}#pi^{-}#pi^{0})", 120, 0.5, 1.0);

        H1F hQ2 = new H1F("hQ2", "Q^{2}", 100, 0.0, 5.0);
        H1F hW  = new H1F("hW", "W", 100, 1.0, 3.0);
        H1F hy  = new H1F("hy", "y", 100, 0.0, 1.0);

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

                int n = rec.getRows();
                if (n < 2) continue;

                ArrayList<Integer> gammas = new ArrayList<>();
                ArrayList<Integer> pips = new ArrayList<>();
                ArrayList<Integer> pims = new ArrayList<>();
                ArrayList<Integer> electrons = new ArrayList<>();

                double vzElectron = 0;

                // =========================
                // PARTICLE SELECTION
                // =========================
                for (int i = 0; i < n; i++) {

                    int pid = rec.getInt("pid", i);
                    int status = rec.getInt("status", i);

                    if (status >= 0) continue;

                    if (pid == 11) electrons.add(i);
                    if (pid == 22) gammas.add(i);
                    if (pid == 211) pips.add(i);
                    if (pid == -211) pims.add(i);
                }

                if (electrons.isEmpty()) continue;
                if (gammas.size() < 2) continue;
                if (pips.isEmpty() || pims.isEmpty()) continue;

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

                vzElectron = rec.getFloat("vz", bestE);

                // =========================
                // VERTEX CUT (TARGET SEPARATION)
                // =========================
                if (vzElectron < -13 || vzElectron > -7) continue;

                // =========================
                // DIS VARIABLES (REAL FORMULAS)
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

                // =========================
                // γγ → π0
                // =========================
                for (int i = 0; i < gammas.size(); i++) {
                    for (int j = i + 1; j < gammas.size(); j++) {

                        LorentzVector g1 = getVector(rec, gammas.get(i));
                        LorentzVector g2 = getVector(rec, gammas.get(j));

                        LorentzVector pi0 = new LorentzVector();
                        pi0.add(g1);
                        pi0.add(g2);

                        double mpi0 = pi0.mass();
                        hPi0.fill(mpi0);

                        // π0 selection
                        if (Math.abs(mpi0 - 0.134) > 0.025) continue;

                        // =========================
                        // ANGLE CUT e- vs γ
                        // =========================
                        for (int g : gammas) {

                            LorentzVector gamma = getVector(rec, g);

                            double angle =
                                    Math.acos(
                                            electron.vect().dot(gamma.vect()) /
                                                    (electron.vect().mag() * gamma.vect().mag())
                                    );

                            if (angle < 0.02) continue; // reject collinear EM brem

                            // =========================
                            // ω → π+ π- π0
                            // =========================
                            for (int ip : pips) {
                                for (int im : pims) {

                                    LorentzVector pip = getVector(rec, ip);
                                    LorentzVector pim = getVector(rec, im);

                                    LorentzVector omega = new LorentzVector();
                                    omega.add(pip);
                                    omega.add(pim);
                                    omega.add(pi0);

                                    double m = omega.mass();

                                    if (m > 0.65 && m < 0.9) {
                                        hOmega.fill(m);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            reader.close();
        }

        // =========================
        // π0 FIT
        // =========================
        F1D fit = new F1D("fit", "[amp]*gaus(x,[mean],[sigma])", 0.10, 0.17);

        fit.setParameter(0, hPi0.getMax());
        fit.setParameter(1, 0.135);
        fit.setParameter(2, 0.01);

        DataFitter.fit(fit, hPi0, "Q");

        System.out.println("\n===== π0 FIT =====");
        System.out.println("Mean  = " + fit.getParameter(1));
        System.out.println("Sigma = " + fit.getParameter(2));

        // =========================
        // PLOTS
        // =========================
        TCanvas c1 = new TCanvas("DIS", 900, 600);
        c1.divide(2,2);
        c1.cd(0); c1.draw(hQ2);
        c1.cd(1); c1.draw(hW);
        c1.cd(2); c1.draw(hy);
        c1.cd(3); c1.draw(hPi0);

        TCanvas c2 = new TCanvas("Omega", 800, 600);
        c2.draw(hOmega);


        System.out.println("Done.");
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