{
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);
    volScalarField rAU2(1.0/UEqn.A());
    volVectorField HbyA2(constrainHbyA(rAU2*UEqn.H(), U, zeroGp));
    surfaceScalarField phiHbyA2("phiHbyA2", fvc::flux(HbyA2));
    MRF.makeRelative(phiHbyA2);
    adjustPhi(phiHbyA2, U, zeroGp);

    tmp<volScalarField> rAtU(rAU2);

    if (simple.consistent())
    {
        rAtU = 1.0/(1.0/rAU2 - UEqn.H1());
        phiHbyA2 +=
            fvc::interpolate(rAtU() - rAU2)*fvc::snGrad(p)*mesh.magSf();
        HbyA2 -= (rAU2 - rAtU())*fvc::grad(zeroGp);
    }

    tUEqn.clear();

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(zeroGp, U, phiHbyA2, rAtU(), MRF);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), zeroGp) == fvc::div(phiHbyA2)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        // if (simple.finalNonOrthogonalIter())
        // {
        //     phi2 = phiHbyA2 - pEqn.flux();
        // }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    // zeroGp.relax();

    // // Momentum corrector
    // U = HbyA2 - rAtU()*fvc::grad(zeroGp);
    // U.correctBoundaryConditions();
    // fvOptions.correct(U);
}
