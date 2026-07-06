import org.jlab.groot.ui.TCanvas;

/**
 * All the TCanvas setup/drawing lives here, separate from the event loop,
 * so OmegaRG_E.main() is just "read events, apply cuts, fill histograms"
 * followed by one call to Plotting.drawAll(h).
 */
public class Plotting {

    public static void drawAll(Histograms h) {

        TCanvas c1 = new TCanvas("DIS", 900, 600);
        c1.divide(2, 2);
        c1.cd(0); c1.draw(h.hQ2);
        c1.cd(1); c1.draw(h.hW);
        c1.cd(2); c1.draw(h.hy);
        c1.cd(3); c1.draw(h.hPi0);

        TCanvas c2 = new TCanvas("Omega", 800, 600);
        c2.draw(h.hOmega);

        TCanvas c3 = new TCanvas("Multiplicities", 900, 600);
        c3.divide(2, 2);
        c3.cd(0); c3.draw(h.hNElec);
        c3.cd(1); c3.draw(h.hNGamma);
        c3.cd(2); c3.draw(h.hNPip);
        c3.cd(3); c3.draw(h.hNPim);

        TCanvas c4 = new TCanvas("Electron_Kinematics", 900, 600);
        c4.divide(2, 2);
        c4.cd(0); c4.draw(h.hElecP);
        c4.cd(1); c4.draw(h.hElecTh);
        c4.cd(2); c4.draw(h.hElecPh);
        c4.cd(3); c4.draw(h.hElecVz);

        TCanvas c5 = new TCanvas("Pip_Kinematics", 900, 600);
        c5.divide(2, 2);
        c5.cd(0); c5.draw(h.hPipP);
        c5.cd(1); c5.draw(h.hPipTh);
        c5.cd(2); c5.draw(h.hPipPh);
        c5.cd(3); c5.draw(h.hPipVz);

        TCanvas c6 = new TCanvas("Pim_Kinematics", 900, 600);
        c6.divide(2, 2);
        c6.cd(0); c6.draw(h.hPimP);
        c6.cd(1); c6.draw(h.hPimTh);
        c6.cd(2); c6.draw(h.hPimPh);
        c6.cd(3); c6.draw(h.hPimVz);

        TCanvas c7 = new TCanvas("Gamma_Kinematics", 900, 600);
        c7.divide(2, 2);
        c7.cd(0); c7.draw(h.hGamP);
        c7.cd(1); c7.draw(h.hGamTh);
        c7.cd(2); c7.draw(h.hGamPh);
        c7.cd(3); c7.draw(h.hGamVz);

        // omega system kinematics + e-gamma angle
        TCanvas c8 = new TCanvas("Omega_Kinematics_and_Angle", 900, 600);
        c8.divide(2, 2);
        c8.cd(0); c8.draw(h.hOmegaP);
        c8.cd(1); c8.draw(h.hOmegaTh);
        c8.cd(2); c8.draw(h.hOmegaPh);
        c8.cd(3); c8.draw(h.hEGammaAngle);
    }
}
