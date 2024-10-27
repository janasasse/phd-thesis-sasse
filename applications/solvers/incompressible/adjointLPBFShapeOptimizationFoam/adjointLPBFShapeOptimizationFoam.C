/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    ajointLPBFShapeOptimizationFoam

Description
    Steady-state solver for incompressible, laminar flow of non-Newtonian
    fluids with optimisation for minimal pressure loss or specified thermal or
    material value at outletusing an adjoint formulation under consideration
    of design restrictions from additive manufacturing (assuming build in +x).

    References:
    \verbatim
        "Adjoint-based topology optimisation of polymer melt flow channels
        producible by additive manufacturing"
        J. Sasse,
        PhD thesis (2024)
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

template<class Type>
void zeroCells
(
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells
)
{
    forAll(cells, i)
    {
        vf[cells[i]] = Zero;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "initAdjointContinuityErrs.H"

    #include "initDataStructure.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop(runTime))
    {

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (initialization != "yes")
        {
            #include "optimiseGeometry.H"
        }

        // Pressure-velocity SIMPLE corrector
        {

            // set local material properties
            forAll(alpha, celli)
            {
                if (alpha[celli] <= (0.5*alphaMax.value()))
                {
                    DT[celli] = DTfluid.value();
                    DConc[celli] = DConcfluid.value();
                    c[celli] = cfluid.value();
                }
                else
                {
                    DT[celli] = DTsolid.value();
                    DConc[celli] = DConcsolid.value();
                    c[celli] = csolid.value();
                }
            }

            // Momentum predictor

            tmp<fvVectorMatrix> tUEqn
            (
                fvm::div(phi, U)
              + turbulence->divDevReff(U)
              + fvm::Sp(alpha, U)
             ==
                fvOptions(U)
            );
            fvVectorMatrix& UEqn = tUEqn.ref();

            UEqn.relax();

            fvOptions.constrain(UEqn);

            solve(UEqn == -fvc::grad(p));

            fvOptions.correct(U);

            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            tUEqn.clear();
            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Momentum corrector
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
            fvOptions.correct(U);

            // Solve primal temperature problem
            // source and documentation: www.dhcae-tools.com/docs/viscousHeatingSolver.pdf
            volTensorField gradU = fvc::grad(U);
            volScalarField nu = laminarTransport.nu();
            tau = nu * (gradU + gradU.T());

            fvScalarMatrix TEqn
            (
                fvm::div(phi, T)
                - fvm::laplacian(DT, T)
                - (1/c)*(tau && gradU)
                ==
                fvOptions(T)
            );

            TEqn.relax();
            TEqn.solve().initialResidual();

            // Solve primal concentration problem
            fvScalarMatrix ConcEqn
            (
                fvm::div(phi,Conc)
                -fvm::laplacian(DConc, Conc)
                ==
                fvOptions(Conc)
            );

            ConcEqn.relax();
            ConcEqn.solve().initialResidual();

        }

        // Adjoint Pressure-velocity SIMPLE corrector
        {
            // Adjoint Momentum predictor

            volVectorField adjointTransposeConvection((fvc::grad(Ua) & U));
            // volVectorField adjointTransposeConvection
            //(
            //    fvc::reconstruct
            //    (
            //        mesh.magSf()*fvc::dotInterpolate(fvc::snGrad(Ua), U)
            //    )
            //);

            zeroCells(adjointTransposeConvection, inletCells);

            volTensorField gradU = fvc::grad(U);
            volVectorField gradTa = fvc::grad(Ta);
            volVectorField gradConca = fvc::grad(Conca);
            volScalarField nu = laminarTransport.nu();
            tau = nu * (gradU + gradU.T());

            tmp<fvVectorMatrix> tUaEqn
            (
                fvm::div(-phi, Ua)
              - adjointTransposeConvection
              + turbulence->divDevReff(Ua)
              + fvm::Sp(alpha, Ua)
              - (T * gradTa)
              + (1/c) * (gradTa & tau)
              - (Conc * gradConca)
             ==
                fvOptions(Ua)
            );
            fvVectorMatrix& UaEqn = tUaEqn.ref();

            UaEqn.relax();

            fvOptions.constrain(UaEqn);

            solve(UaEqn == -fvc::grad(pa));

            fvOptions.correct(Ua);

            volScalarField rAUa(1.0/UaEqn.A());
            volVectorField HbyAa("HbyAa", Ua);
            HbyAa = rAUa*UaEqn.H();
            tUaEqn.clear();
            surfaceScalarField phiHbyAa("phiHbyAa", fvc::flux(HbyAa));
            adjustPhi(phiHbyAa, Ua, pa);

            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix paEqn
                (
                    fvm::laplacian(rAUa, pa) == fvc::div(phiHbyAa)
                );

                paEqn.setReference(paRefCell, paRefValue);
                paEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phia = phiHbyAa - paEqn.flux();
                }
            }

            #include "adjointContinuityErrs.H"

            // Explicitly relax pressure for adjoint momentum corrector
            pa.relax();

            // Adjoint momentum corrector
            Ua = HbyAa - rAUa*fvc::grad(pa);
            Ua.correctBoundaryConditions();
            fvOptions.correct(Ua);

            // solve dual problem for temperature
            fvScalarMatrix TaEqn
            (
                fvm::div(-phi, Ta)
                - fvm::laplacian(DT, Ta)
                ==
                fvOptions(Ta)
            );

            TaEqn.relax();
            fvOptions.constrain(TaEqn);
            TaEqn.solve().initialResidual();
            fvOptions.correct(Ta);

            // Solve dual concentration problem
            forAll (alpha, celli)
            {
                DConc[celli] = DConca.value();
            }

            fvScalarMatrix ConcaEqn
            (
                fvm::div(-phi, Conca)
                - fvm::laplacian(DConc, Conca)
                ==
                fvOptions(Conca)
            );

            ConcaEqn.relax();
            ConcaEqn.solve().initialResidual();

        }

        laminarTransport.correct();
        turbulence->correct();

        sens = (
            lambdaU*(Ua & U)
          + lambdaT*wT*(Ta * (T - average(T)))
          + lambdaC*wC*(Conca * (Conc - average(Conc))));

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
