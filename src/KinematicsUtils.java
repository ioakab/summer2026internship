import org.jlab.clas.physics.LorentzVector;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.groot.data.H1F;

/**
 * Small collection of static helper methods shared by OmegaRG_E:
 *  - building a LorentzVector from a REC::Particle row
 *  - computing the opening angle between two particles
 *  - filling the standard p/theta/phi/vz histogram set for one particle
 *  - proton mass constant
 */
public class KinematicsUtils {

    // proton mass
    public static double M_p() {
        return 0.938272;
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

    // =========================
    // Helper: opening angle between two particles (radians)
    // =========================
    public static double angleBetween(LorentzVector a, LorentzVector b) {
        return Math.acos(
                a.vect().dot(b.vect()) / (a.vect().mag() * b.vect().mag())
        );
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
}
