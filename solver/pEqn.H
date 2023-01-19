/*Switch CombinedHSPressure("no");
runTime.controlDict().readIfPresent("CombinedHSPressure", CombinedHSPressure);
const scalar CHSP = CombinedHSPressure ? 1.0 : 0;

{
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));

    U = rAU()*UEqn.H();


    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
        + fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi)
    );
   MRF.makeRelative(phiU);

    adjustPhi(phiU, U, p_rgh);

    //Calculate the flux from surface tension:
    surfaceScalarField phi_c
    (
        IOobject
        (
            "phi_c",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("flux", dimLength*dimLength*dimLength/dimTime, 0)
    );


    phi_c   = stfModel->phi_c( rAUf );      

    //Calculate the flux from gravitational forces
    surfaceScalarField phi_rgh
    (
        IOobject
        (
            "phi_rgh",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("flux", dimLength*dimLength*dimLength/dimTime, 0)
    );

    if (!CombinedHSPressure)
    {
        phi_rgh =  -ghf*fvc::snGrad(rho)*rAUf*mesh.magSf();
    }

    phi = phiU + stfModel->phi_c( rAUf ) + phi_rgh;


   //surfaceScalarField phiHbyA = phi_rgh.boundaryField() - MRF.relative(mesh.Sf().boundaryField() & U.boundaryField()))/(mesh.magSf().boundaryField()*rAUf.boundaryField())

   // Update the fixedFluxPressure BCs to ensure flux consistency
    //constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);
    
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryFieldRef(),
        (
            phi_rgh.boundaryField() - MRF.relative(mesh.Sf().boundaryField() & U.boundaryField())
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    while (pimple.correctNonOrthogonal())
    {

    // No need to correct p_rgh for phase change
    // - already incorporated into div(phi) 
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh)
            == 
            fvc::div(phi)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

    // This is -= (unlike phaseChangeFoam) because the sign of the governing
    // equations is changed.
        if (pimple.finalNonOrthogonalIter())
        {
            phi -= p_rghEqn.flux();
        }

    }
    Info<< "\nTest1\n" << endl;
    // Added Condition for Modified Surface Tension Case 
    U += rAU()*fvc::reconstruct((phi - phiU)/rAUf); 

    //U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
    Info<< "\nTest2\n" << endl;
    U.correctBoundaryConditions();

    #include "continuityErrs.H"
    //Add in surface tension
    p == p_rgh + (1.0-CHSP)*rho*gh + stfModel->pcap();
    Info<< "\nTest3\n" << endl;
    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
           "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        //Add in surface tension
        p_rgh = p - (1.0-CHSP)*rho*gh - stfModel->pcap();
    }

}

*/



{
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
    );
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    if (!correctPhi)
    {
        rAU.clear();
    }
}
