import org.jlab.clas.physics.LorentzVector;
import org.jlab.jnp.hipo4.io.*;
import org.jlab.jnp.hipo4.data.*;
import java.util.ArrayList;

void main() {
    // ----------------------------
    // 1. Open file
    // ----------------------------
    HipoReader reader = new HipoReader();
    reader.open("input.hipo"); // <-- change this

    Event event = new Event();
    SchemaFactory factory = reader.getSchemaFactory();
    Bank particles = new Bank(factory.getSchema("REC::Particle"));

    // simple histogram bins (ω mass)
    double[] hist = new double[200];
    double min = 0.5;
    double max = 1.2;

    // ----------------------------
    // 2. Event loop
    // ----------------------------
    while (reader.hasNext()) {

        reader.nextEvent(event);
        event.read(particles);

        int n = particles.getRows();
        if (n == 0) continue;

        // ----------------------------
        // 3. Store particles
        // ----------------------------
        ArrayList<int[]> gammas = new ArrayList<>();
        ArrayList<int[]> pips = new ArrayList<>();
        ArrayList<int[]> pims = new ArrayList<>();

        for (int i = 0; i < n; i++) {

            int pid = particles.getInt("pid", i);

            if (pid == 22) {
                gammas.add(new int[]{i});
            }
            if (pid == 211) {
                pips.add(new int[]{i});
            }
            if (pid == -211) {
                pims.add(new int[]{i});
            }
        }

        // need at least 2 gammas + π+ + π-
        if (gammas.size() < 2 || pips.size() < 1 || pims.size() < 1)
            continue;

        // ----------------------------
        // 4. Build π0 from γγ pairs
        // ----------------------------
        for (int i = 0; i < gammas.size(); i++) {
            for (int j = i + 1; j < gammas.size(); j++) {

                LorentzVector g1 = getVector(particles, gammas.get(i)[0]);
                LorentzVector g2 = getVector(particles, gammas.get(j)[0]);

                LorentzVector pi0 = new LorentzVector();
                pi0.add(g1);
                pi0.add(g2);

                // optional π0 mass cut
                if (pi0.mass() < 0.08 || pi0.mass() > 0.20)
                    continue;

                // ----------------------------
                // 5. Build ω = π+ π- π0
                // ----------------------------
                for (int ip = 0; ip < pips.size(); ip++) {
                    for (int im = 0; im < pims.size(); im++) {

                        LorentzVector pip = getVector(particles, pips.get(ip)[0]);
                        LorentzVector pim = getVector(particles, pims.get(im)[0]);

                        LorentzVector omega = new LorentzVector();
                        omega.add(pip);
                        omega.add(pim);
                        omega.add(pi0);

                        double m = omega.mass();

                        // fill histogram
                        int bin = (int) ((m - min) / (max - min) * hist.length);
                        if (bin >= 0 && bin < hist.length)
                            hist[bin]++;
                    }
                }
            }
        }
    }

    reader.close();

    // ----------------------------
    // 6. Print histogram
    // ----------------------------
    for (int i = 0; i < hist.length; i++) {
        double x = min + i * (max - min) / hist.length;
        System.out.println(x + " " + hist[i]);
    }
}

// ----------------------------
// helper: build 4-vector
// ----------------------------
public static LorentzVector getVector(Bank b, int row) {

    double px = b.getFloat("px", row);
    double py = b.getFloat("py", row);
    double pz = b.getFloat("pz", row);
    double e  = b.getFloat("energy", row);

    LorentzVector v = new LorentzVector();
    v.setPxPyPzE(px, py, pz, e);

    return v;
}

