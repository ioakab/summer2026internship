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

public class Main {

    public static void main(String[] args) {

        String path = "/home/teo/Documents/HipoFold";

        File folder = new File(path);
        File[] files = folder.listFiles((d, name) -> name.endsWith(".hipo"));

        if (files == null || files.length == 0) {
            System.out.println("No files found.");
            return;
        }

        // =====================
        // HISTOGRAMS
        // =====================
        H1F hPi0 = new H1F("hPi0", "M(gamma gamma)", 120, 0.0, 0.30);
        H1F hOmega = new H1F("hOmega", "M(pi+ pi- pi0)", 120, 0.5, 1.0);

        // =====================
        // LOOP FILES
        // =====================
        for (File file : files) {

            System.out.println("Processing: " + file.getName());

            HipoReader reader = new HipoReader();
            reader.open(file.getAbsolutePath());

            Event event = new Event();
            SchemaFactory factory = reader.getSchemaFactory();
            Bank rec = new Bank(factory.getSchema("REC::Particle"));

            // =====================
            // EVENT LOOP
            // =====================
            while (reader.hasNext()) {

                reader.nextEvent(event);
                event.read(rec);

                int n = rec.getRows();
                if (n < 2) continue;

                ArrayList<Integer> gammas = new ArrayList<>();
                ArrayList<Integer> pips = new ArrayList<>();
                ArrayList<Integer> pims = new ArrayList<>();
                ArrayList<Integer> electrons = new ArrayList<>();

                // =====================
                // PARTICLE SELECTION
                // =====================
                for (int i = 0; i < n; i++) {

                    int pid = rec.getInt("pid", i);
                    int status = rec.getInt("status", i);

                    if (status >= 0) continue;

                    if (pid == 22) gammas.add(i);
                    else if (pid == 211) pips.add(i);
                    else if (pid == -211) pims.add(i);
                    else if (pid == 11) electrons.add(i);
                }

                if (gammas.size() < 2) continue;
                if (pips.isEmpty() || pims.isEmpty()) continue;
                if (electrons.isEmpty()) continue;

                // =====================
                // DIS CUT (placeholder)
                // =====================
                LorentzVector e = getVector(rec, electrons.get(0));
                double Q2 = e.px()*e.px() + e.py()*e.py() + e.pz()*e.pz();

                if (Q2 < 1.0) continue;

                // =====================
                // π0 reconstruction
                // =====================
                for (int i = 0; i < gammas.size(); i++) {
                    for (int j = i + 1; j < gammas.size(); j++) {

                        LorentzVector g1 = getVector(rec, gammas.get(i));
                        LorentzVector g2 = getVector(rec, gammas.get(j));

                        LorentzVector pi0 = new LorentzVector();
                        pi0.add(g1);
                        pi0.add(g2);

                        double mPi0 = pi0.mass();
                        hPi0.fill(mPi0);

                        // π0 cut
                        if (Math.abs(mPi0 - 0.134) > 0.025) continue;

                        // =====================
                        // ω reconstruction
                        // =====================
                        for (int ip : pips) {
                            for (int im : pims) {

                                LorentzVector pip = getVector(rec, ip);
                                LorentzVector pim = getVector(rec, im);

                                LorentzVector omega = new LorentzVector();
                                omega.add(pip);
                                omega.add(pim);
                                omega.add(pi0);

                                double mOmega = omega.mass();

                                if (mOmega > 0.65 && mOmega < 0.9) {
                                    hOmega.fill(mOmega);
                                }
                            }
                        }
                    }
                }
            }

            reader.close();
        }

        // =====================
        // π0 FIT
        // =====================
        F1D fit = new F1D("fit", "[amp]*gaus(x,[mean],[sigma])", 0.10, 0.17);

        fit.setParameter(0, hPi0.getMax());
        fit.setParameter(1, 0.135);
        fit.setParameter(2, 0.01);

        DataFitter.fit(fit, hPi0, "Q");

        System.out.println("Pi0 mean: " + fit.getParameter(1));
        System.out.println("Pi0 sigma: " + fit.getParameter(2));

        // =====================
        // DRAW
        // =====================
        TCanvas c1 = new TCanvas("pi0", 800, 600);
        c1.draw(hPi0);
        c1.draw(fit, "same");

        TCanvas c2 = new TCanvas("omega", 800, 600);
        c2.draw(hOmega);

        System.out.println("Done.");
    }

    // =====================
    // LORENTZ VECTOR BUILDER
    // =====================
    public static LorentzVector getVector(Bank b, int row) {

        double px = b.getFloat("px", row);
        double py = b.getFloat("py", row);
        double pz = b.getFloat("pz", row);

        int pid = b.getInt("pid", row);

        double mass;

        if (pid == 22) mass = 0.0;
        else if (pid == 211 || pid == -211) mass = 0.13957;
        else if (pid == 11) mass = 0.000511;
        else mass = 0.0;

        LorentzVector v = new LorentzVector();
        v.setPxPyPzM(px, py, pz, mass);

        return v;
    }
}