/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | avalanche module
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Matthias Rauter
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
    faParkerFukushimaFoam

Description
    A depth-integrated solver for turbidity currents.
    The solver is based on the Finite Area Method.
    Model derivation and description:
    Parker,, Fukushima, and Pantin (1986). "Self-accelerating
    turbidity currents". Journal of Fluid Mechanics, 171, 145-181.
    dx.doi.org/10.1017/S0022112086001404

Author
    Matthias Rauter matthias@rauter.it

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"
#include "SolverPerformance.H"
#include "suspensionFrictionModel.H"
#include "suspensionEntrainmentModel.H"
#include "suspensionDepositionModel.H"
#include "ambientEntrainmentModel.H"

#define CREATE_FIELDS createFaFields.H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "(avalanche)\n"
        "A depth-integrated solver for turbidity currents using"
        " Finite Area Methods."
    );

    #include "postProcessFA.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFaMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "readTransportProperties.H"
    #include "createFaFields.H"
    #include "createSubmodels.H"
    #include "createTimeControls.H"


    Info<< "\nStarting time loop\n" << endl;

    #include "readSolutionControls.H"

    Info<< nl
        << "Numerical settings" << nl
        << "    max number of iterations " << nCorr << nl
        << "    min number of iterations " << minCorr << nl
        << "    TOL Us " << UsResidualMax << nl
        << "    TOL h " << hResidualMax << nl
        << "    TOL c " << cResidualMax << nl << endl;

    const bool initDeltaT = runTime.controlDict().get<bool>("initDeltaT");

    if (initDeltaT)
    {
        Info<< "Initializing Delta T" << endl;
        #include "readTimeControls.H"
        #include "surfaceCourantNo.H"
        runTime.setDeltaT
        (
            min(maxCo/(CoNum + SMALL)*runTime.deltaT().value(), maxDeltaT)
        );
    }

    bool final = false;
    bool success = false;

    while (runTime.run())
    {
        #include "readSolutionControls.H"
        #include "readTimeControls.H"
        #include "surfaceCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        final = false;

        for (int iCorr = 0; ; iCorr++)
        {

            Info<< "Correction iteration " << iCorr << endl;

            geff = max(gn - (fac::ndiv(phis, Us)&n), dimensionedScalar(dimVelocity/dimTime));

            const areaVectorField & tauSc = friction->tauSc();
            const areaScalarField & tauSp = friction->tauSp();

            //Parker et al. (1986), Eq. (5)
            faVectorMatrix UsEqn
            (
                fam::ddt(h, Us)
              + fam::div(phi2s, Us)
              + tauSc               // ustar**2 expl. part
              + fam::Sp(tauSp, Us)  // ustar**2 impl. part
              ==
                R*gs*c*h
              - 1./2.*R*geff*fac::grad(c*h*h)*boundaryCell
            );

            if (correctMomentum)
            {
                UsEqn += R*c*fam::ddt(h, Us)
                       + R*h*fam::ddt(c, Us)
                       + R*fam::div(fac::interpolate(c)*phi2s, Us);
            }

            if (!final)
                UsEqn.relax();

            SolverPerformance<vector> UsResidual = solve(UsEqn);

            Us.correctBoundaryConditions();

            phis = (fac::interpolate(Us) & aMesh.Le());

            tau = tauSc + tauSp*Us;



            //Calculate water entrainment
            Sa = ambientEntrainment->Sm();
            Sa = Sa*pos(hwem-h); //clip entrainment if maximum water height is reached

            //Additional sink term for water
            const areaScalarField vs(R*geff*Ds*Ds/18./nu); //Terminal velocity
            areaScalarField dw(pos(h-hmin*2.)*vs);


            //Calculate sediment entrainment
            areaScalarField Sm = entrainment->Sm();
            //no entrainment if there is not a minimal flow thickness
            Sm = Sm*pos(h-hentmin);
            Sm = min(Sm, hentrain/runTime.deltaT());


            //Parker et al. (1986), Eq. (3)
            faScalarMatrix hEqn
            (
                fam::ddt(h)
              + fam::div(phis, h)
              ==
                Sa
//              + Sm/(hentrain+dimensionedScalar("hentrainmin", dimLength, 1e-5))*hentrain
            );

            if (waterSink)
            {
                hEqn += fam::Sp(dw/(h+dimensionedScalar("hmin", dimLength, 1e-5)),h);
            }

            solverPerformance hResidual = hEqn.solve();

            phi2s = hEqn.flux();

            if (!final)
            {
                 h.relax();
            }

            if (bindHeight)
            {
                h = max(h, hmin);
            }

            //Calculate sedimation rate
            Sd = sedimentation->Sd();
            Sd = min(Sd, h*c/runTime.deltaT());


            faScalarMatrix cEqn
            (
                fam::ddt(h, c)
              + fam::div(phi2s, c)
             ==
                Sm/(hentrain+dimensionedScalar("hentrainmin", dimLength, 1e-5))*hentrain
              - fam::Sp(Sd/(c+dimensionedScalar("cmin", dimless, 1e-8)), c)
            );

            if(!final)
                cEqn.relax();

            solverPerformance cResidual = cEqn.solve();

            if (bindHeight)
                c = max(c, cmin);

            //Keep track of bed sediment
            faScalarMatrix hentrainEqn
            {
                fam::ddt(hentrain)
             ==
              - fam::Sp(Sm/(hentrain+dimensionedScalar("heintrainmin", dimLength, 1e-5)), hentrain)
              + Sd/(c+dimensionedScalar("cmin", dimless, 1e-8))*c
            };

            hentrainEqn.solve();


            if (final)
            {
                if
                (
                    success
                 && mag(UsResidual.initialResidual()) < UsResidualMax
                 && hResidual.initialResidual() < hResidualMax
                 && cResidual.initialResidual() < cResidualMax
                )
                {
                    Info<< "Reached residual in Us = "
                        << UsResidual.initialResidual()
                        << " < " << UsResidualMax
                        << ", in h = "
                        << hResidual.initialResidual()
                        << " < " << hResidualMax
                        << ", in c = "
                        << cResidual.initialResidual()
                        << " < " << cResidualMax
                        << ", stopping loop!" << endl;
                     break;
                }
                if (!success)
                {
                    Info<< "Reached maximum number of iterations (" << iCorr << "). Residual in Us = "
                        << UsResidual.initialResidual()
                        << ", in h = "
                        << hResidual.initialResidual()
                        << ", in c = "
                        << cResidual.initialResidual()
                        << ", stopping loop!" << endl;
                    break;
                }
            }

            if
            (
                (
                    mag(UsResidual.initialResidual()) < UsResidualMax
                 && hResidual.initialResidual() < hResidualMax
                 && cResidual.initialResidual() < cResidualMax
                 && iCorr >= minCorr
                )
             )
            {
                final = true;
                success = true;
            }
            else if (iCorr >= nCorr)
            {
                final = true;
                success = false;
            }
        }

        if (runTime.outputTime())
        {
            runTime.write();
        }

        runTime.printExecutionTime(Info);
    }


    Info<< nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
