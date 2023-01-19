/*
// Combined pressure model
Switch CombinedHSPressure("no");
runTime.controlDict().readIfPresent("CombinedHSPressure", CombinedHSPressure);
const scalar CHSP = CombinedHSPressure ? 1.0 : 0;

// Improved transport models
Switch ImprovedTransportBlending("no");
transportProperties.readIfPresent("ImprovedTransportBlending", ImprovedTransportBlending);
const scalar ImpTB = ImprovedTransportBlending ? 1.0 : 0.0;    
*/
MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
      //ASR - the default divDevRhoReff uses the viscosity values on cell centers
      // better to use the face values with this improved blending model
  //    + ImpTB*(   -fvm::laplacian(mixture.muf(), U)
  //             + fvm::laplacian(mixture.mu(), U) )
  //    - CHSP*rho*g
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
     //               stfModel->Fstff()
                   mixture.surfaceTensionForce()
     //             - (1.0-CHSP)*ghf*fvc::snGrad(rho)
                    - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
