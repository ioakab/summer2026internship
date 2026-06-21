import org.jlab.clas.physics.LorentzVector;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.SchemaFactory;

import org.jlab.groot.data.H1F;
import org.jlab.groot.ui.TCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.fitter.DataFitter;

import java.util.ArrayList;
import java.util.Arrays;
import java.io.File;

void main() {

    // =================================================
    // HIPO directory
    // =================================================
    String directoryPath = "/home/teo/Documents/HipoFold"; // CHANGE THIS

    File folder = new File(directoryPath);

    File[] hipoFiles = folder.listFiles(
            (dir, name) -> name.endsWith(".hipo"));

    if (hipoFiles == null || hipoFiles.length == 0) {
        System.out.println("No HIPO files found.");
        return;
    }

    Arrays.sort(hipoFiles);

    // =================================================
    // Histograms
    // =================================================
    H1F hPi0 = new H1F("hPi0",
            "Invariant Mass #gamma#gamma",
            120, 0.0, 0.30);

    hPi0.setTitleX("M(#gamma#gamma) [GeV]");
    hPi0.setTitleY("Counts");

    H1F hOmega = new H1F("hOmega",
            "#omega → π⁺ π⁻ π⁰",
            120, 0.5, 1.0);

    hOmega.setTitleX("M(π⁺π⁻π⁰) [GeV]");
    hOmega.setTitleY("Counts");

    // =================================================
    // Loop over files
    // =================================================
    for (File file : hipoFiles) {

        System.out.println("Processing: " + file.getName());

        HipoReader reader = new HipoReader();
        reader.open(file.getAbsolutePath());

        Event event = new Event();
        SchemaFactory factory = reader.getSchemaFactory();
        Bank particles =
                new Bank(factory.getSchema("REC::Particle"));

        // =============================================
        // Event loop
        // =============================================
        while (reader.hasNext()) {

            reader.nextEvent(event);
            event.read(particles);

            int n = particles.getRows();
            if (n == 0) continue;

            ArrayList<Integer> gammas = new ArrayList<>();
            ArrayList<Integer> pips = new ArrayList<>();
            ArrayList<Integer> pims = new ArrayList<>();
            ArrayList<Integer> electrons = new ArrayList<>();

            // =========================================
            // STEP 1 (CRITICAL FIX): CLAS12 status cut
            // =========================================
            for (int i = 0; i < n; i++) {

                int pid = particles.getInt("pid", i);
                int status = particles.getInt("status", i);

                if (status >= 0) continue; // MUST HAVE THIS

                if (pid == 22) gammas.add(i);
                if (pid == 211) pips.add(i);
                if (pid == -211) pims.add(i);
                if (pid == 11) electrons.add(i);
            }

            // Require final state
            if (gammas.size() < 2) continue;
            if (pips.size() < 1) continue;
            if (pims.size() < 1) continue;
            if (electrons.size() < 1) continue;

            // =========================================
            // γγ → π0
            // =========================================
            for (int i = 0; i < gammas.size(); i++) {
                for (int j = i + 1; j < gammas.size(); j++) {

                    LorentzVector g1 = getVector(particles, gammas.get(i));
                    LorentzVector g2 = getVector(particles, gammas.get(j));

                    LorentzVector pi0 = new LorentzVector();
                    pi0.add(g1);
                    pi0.add(g2);

                    double mpi0 = pi0.mass();

                    hPi0.fill(mpi0);

                    // π0 mass cut
                    if (Math.abs(mpi0 - 0.134) > 0.025)
                        continue;

                    // =====================================
                    // ω → π⁺ π⁻ π⁰
                    // =====================================
                    for (int ip = 0; ip < pips.size(); ip++) {
                        for (int im = 0; im < pims.size(); im++) {

                            LorentzVector pip = getVector(particles, pips.get(ip));
                            LorentzVector pim = getVector(particles, pims.get(im));

                            LorentzVector omega = new LorentzVector();
                            omega.add(pip);
                            omega.add(pim);
                            omega.add(pi0);

                            double m = omega.mass();

                            // optional ω window (VERY IMPORTANT)
                            if (m > 0.65 && m < 0.9)
                                hOmega.fill(m);
                        }
                    }
                }
            }
        }

        reader.close();
    }

    // =================================================
    // π0 fit
    // =================================================
    F1D fPi0 = new F1D(
            "fPi0",
            "[amp]*gaus(x,[mean],[sigma])",
            0.10, 0.17);

    fPi0.setParameter(0, hPi0.getMax());
    fPi0.setParameter(1, 0.135);
    fPi0.setParameter(2, 0.01);

    DataFitter.fit(fPi0, hPi0, "Q");

    System.out.println("\n===== PI0 FIT =====");
    System.out.println("Mean  = " + fPi0.getParameter(1));
    System.out.println("Sigma = " + fPi0.getParameter(2));

    // =================================================
    // Draw
    // =================================================
    TCanvas c1 = new TCanvas("Pi0", 800, 600);
    c1.draw(hPi0);
    c1.draw(fPi0, "same");

    TCanvas c2 = new TCanvas("Omega", 800, 600);
    c2.draw(hOmega);

    System.out.println("Done.");
}

// =====================================================
// Lorentz vector builder (FIXED for CLAS12)
// =====================================================
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