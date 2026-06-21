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

void main() {

    // -------------------------------------------------
    // Open HIPO file
    // -------------------------------------------------
    HipoReader reader = new HipoReader();
    reader.open("input.hipo");   // Replace with your file name

    Event event = new Event();
    SchemaFactory factory = reader.getSchemaFactory();
    Bank particles = new Bank(factory.getSchema("REC::Particle"));

    // -------------------------------------------------
    // Histograms
    // -------------------------------------------------

    // Gamma-Gamma invariant mass (pi0)
    H1F hPi0 = new H1F("hPi0",
            "Invariant Mass of Two Photons",
            100, 0.0, 0.30);

    hPi0.setTitleX("M(#gamma#gamma) [GeV]");
    hPi0.setTitleY("Counts");

    // Omega invariant mass
    H1F hOmega = new H1F("hOmega",
            "#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}",
            100, 0.50, 1.10);

    hOmega.setTitleX("M(#pi^{+}#pi^{-}#pi^{0}) [GeV]");
    hOmega.setTitleY("Counts");

    // -------------------------------------------------
    // Event Loop
    // -------------------------------------------------
    while (reader.hasNext()) {

        reader.nextEvent(event);
        event.read(particles);

        int n = particles.getRows();
        if (n == 0) continue;

        // ---------------------------------------------
        // Particle lists
        // ---------------------------------------------
        ArrayList<Integer> gammas = new ArrayList<>();
        ArrayList<Integer> pips = new ArrayList<>();
        ArrayList<Integer> pims = new ArrayList<>();
        ArrayList<Integer> electrons = new ArrayList<>();

        for (int i = 0; i < n; i++) {

            int pid = particles.getInt("pid", i);

            if (pid == 22) gammas.add(i);      // photon
            if (pid == 211) pips.add(i);       // pi+
            if (pid == -211) pims.add(i);      // pi-
            if (pid == 11) electrons.add(i);   // electron
        }

        // ------------------------------------------------
        // Require final state:
        // e + pi+ + pi- + 2 gamma
        // ------------------------------------------------
        if (gammas.size() < 2) continue;
        if (pips.size() < 1) continue;
        if (pims.size() < 1) continue;
        if (electrons.size() < 1) continue;

        // ------------------------------------------------
        // Loop over all gamma-gamma combinations
        // ------------------------------------------------
        for (int i = 0; i < gammas.size(); i++) {

            for (int j = i + 1; j < gammas.size(); j++) {

                LorentzVector g1 =
                        getVector(particles, gammas.get(i));

                LorentzVector g2 =
                        getVector(particles, gammas.get(j));

                // Construct pi0 candidate
                LorentzVector pi0 = new LorentzVector();
                pi0.add(g1);
                pi0.add(g2);

                double mpi0 = pi0.mass();

                // Fill gamma-gamma histogram
                hPi0.fill(mpi0);

                // ------------------------------------------------
                // Pi0 mass cut
                // ------------------------------------------------
                if (mpi0 < 0.115 || mpi0 > 0.155)
                    continue;

                // ------------------------------------------------
                // Build omega from pi+, pi-, pi0
                // ------------------------------------------------
                for (int ip = 0; ip < pips.size(); ip++) {

                    for (int im = 0; im < pims.size(); im++) {

                        LorentzVector pip =
                                getVector(particles, pips.get(ip));

                        LorentzVector pim =
                                getVector(particles, pims.get(im));

                        LorentzVector omega =
                                new LorentzVector();

                        omega.add(pip);
                        omega.add(pim);
                        omega.add(pi0);

                        hOmega.fill(omega.mass());
                    }
                }
            }
        }
    }

    reader.close();

    // -------------------------------------------------
    // Fit pi0 histogram
    // -------------------------------------------------
    F1D fPi0 = new F1D(
            "fPi0",
            "[amp]*gaus(x,[mean],[sigma])",
            0.10, 0.17);

    fPi0.setParameter(0, 1000);
    fPi0.setParameter(1, 0.135);
    fPi0.setParameter(2, 0.01);

    DataFitter.fit(fPi0, hPi0, "Q");

    System.out.println("\n===== PI0 FIT RESULTS =====");
    System.out.println("Mean  = " +
            fPi0.getParameter(1) + " GeV");

    System.out.println("Sigma = " +
            fPi0.getParameter(2) + " GeV");

    // -------------------------------------------------
    // Draw Histograms
    // -------------------------------------------------
    TCanvas c1 = new TCanvas("Pi0", 800, 600);
    c1.draw(hPi0);
    c1.draw(fPi0, "same");

    TCanvas c2 = new TCanvas("Omega", 800, 600);
    c2.draw(hOmega);
}


// =====================================================
// Helper function: build Lorentz vector
// =====================================================
public static LorentzVector getVector(Bank b, int row) {

    double px = b.getFloat("px", row);
    double py = b.getFloat("py", row);
    double pz = b.getFloat("pz", row);
    double e  = b.getFloat("energy", row);

    LorentzVector v = new LorentzVector();
    v.setPxPyPzE(px, py, pz, e);

    return v;
}