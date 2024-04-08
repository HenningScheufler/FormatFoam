/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | avalanche module
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Matthias Rauter
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
    faSavageHutterFoam

Description
    A depth-integrated solver for mixed snow avalanches.
    The dense flow is described with the Savage-Hutter model.
    The powder cloud is described with the Parker-Fukushima model.
    The coupling model is new.

    Publication in preparation.

Author
    Matthias Rauter matthias@rauter.it

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"
#include "frictionModel.H"
#include "entrainmentModel.H"
#include "depositionModel.H"
#include "suspensionFrictionModel.H"
#include "suspensionEntrainmentModel.H"
#include "ambientEntrainmentModel.H"
#include "suspensionDepositionModel.H"
#include "couplingModel.H"
#include "SolverPerformance.H"

#define CREATE_FIELDS SHcreateFaFields.H
#define CREATE_FIELDS_2 PFcreateFaFields.H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "(avalanche)\n"
        "A depth-integrated solver for mixed snow avalanches using"
        " Finite Area Methods."
    );

    #include "postProcessFA.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFaMesh.H"
    #include "readGravitationalAcceleration.H"

    #include "SHreadTransportProperties.H"
    #include "PFreadTransportProperties.H"

    #include "SHcreateFaFields.H"
    #include "PFcreateFaFields.H"

    #include "SHcreateSubmodels.H"
    #include "PFcreateSubmodels.H"

    #include "createTimeControls.H"

    Info<< "\nStarting time loop\n" << endl;

    #include "readSolutionControls.H"

    Info<< nl
        << "Numerical settings" << nl
        << "    max number of iterations " << nCorr << nl
        << "    min number of iterations " << minCorr << nl
        << "    TOL h " << hResidualMax << nl
        << "    TOL Us " << UsResidualMax << nl 
        << "    TOL c " << cResidualMax << nl << endl;

    const bool initDeltaT = runTime.controlDict().get<bool>("initDeltaT");

    if (initDeltaT)
    {
        Info<< "Initializing Delta T" << endl;
        #include "readTimeControls.H"
        #include "surfaceCourantNo1.H"
        #include "surfaceCourantNo2.H"

        runTime.setDeltaT
        (
            min(maxCo/(max(CoNum1, CoNum2) + SMALL)*runTime.deltaT().value(), maxDeltaT)
        );
    }

    bool final = false;

    while (runTime.run())
    {
        #include "readSolutionControls.H"
        #include "readTimeControls.H"
        #include "surfaceCourantNo1.H"
        #include "surfaceCourantNo2.H"

        scalar CoNum = max(CoNum1, CoNum2);
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        final = false;

        for (int iCorr = 0; ; iCorr++)
        {

            //Calculate coupling fluxes
            Sdp = coupling->Sdp();
            Sdp = min(Sdp, h1/runTime.deltaT());

            //Calculate sedimation rate
            Spd = sedimentation2->Sd();
            Spd = min(Spd, h2*c2/runTime.deltaT());

            //Info << "Upward mass flux = [" << min(Sdp).value() << ", " << max(Sdp).value() << "]" << endl;
            //Info << "Upward momentum flux = [" << min(Sdp*mag(Us1)).value() << ", " << max(Sdp*mag(Us1)).value() << "]" << endl;

            powderLayerAlone = pos(hmin1*3-h1);

            #include "advanceDenseFlowModel.H"

            #include "advancePowderCloudModel.H"

            //Keep track of intact snow cover
            faScalarMatrix hentrainEqn
            (
                fam::ddt(hentrain)
              ==
                Sd1
                - fam::Sp
                (
                    Sm1/(hentrain + dimensionedScalar("small", dimLength, SMALL)),
                    hentrain
                )
              - fam::Sp(
                    Sm2/(hentrain+dimensionedScalar("small", dimLength, SMALL))/packingd,
                    hentrain)

            );

            hentrainEqn.solve();



            if (final)
            {
                if
                (
                    h1Residual.initialResidual() < hResidualMax
                 && mag(Us1Residual.initialResidual()) < UsResidualMax
                 && h2Residual.initialResidual() < hResidualMax
                 && mag(Us2Residual.initialResidual()) < UsResidualMax
                 && c2Residual.initialResidual() < cResidualMax
                )
                {
                    Info<< "Layer 1: Reached residual in h1 = "
                        << h1Residual.initialResidual()
                        << " < " << hResidualMax
                        << " and in Us1 = "
                        << Us1Residual.initialResidual()
                        << " < " << UsResidualMax << endl
                        << "Layer 2: Reached residual in h2 = "
                        << h2Residual.initialResidual()
                        << " < " << hResidualMax
                        << " and in Us2 = "
                        << Us2Residual.initialResidual()
                        << " < " << UsResidualMax
                        << ", in c2 = "
                        << c2Residual.initialResidual()
                        << " < " << cResidualMax << ", " << endl
                        << "stopping loop!" << endl;
                }
                else
                {
                    Info<< "Reached maximum numbers of iterations, "
                        << "stopping loop!" << endl;
                }

                const scalar Vsl = gSum((hentrain*aMesh.S())())/packingd.value();
                const scalar Vdl = gSum((h1*aMesh.S())())/packingd.value();
                const scalar Vpl = gSum((h2*c2*aMesh.S())());
                Info << "Grain vol. in staionary layer = " << Vsl << endl;
                Info << "Grain vol. in dense flow layer = " << Vdl << endl;
                Info << "Grain vol. in powder cloud layer = " << Vpl << endl;
                Info << "Grain vol. in all layers = " << Vsl + Vdl + Vpl << endl;
                break;
            }

            if
            (
                (
                    h1Residual.initialResidual() < hResidualMax
                 && mag(Us1Residual.initialResidual()) < UsResidualMax
                 && iCorr >= minCorr
                )
             || iCorr >= nCorr
            )
            {
                final = true;
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
