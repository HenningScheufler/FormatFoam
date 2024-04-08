/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | avalanche module
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 Matthias Rauter
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

Author
    Matthias Rauter matthias@rauter.it

\*---------------------------------------------------------------------------*/

#include "suspensionDepositionModel.H"
#include "fvCFD.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(suspensionDepositionModel, 0);
    defineRunTimeSelectionTable(suspensionDepositionModel, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::suspensionDepositionModel::readDict
(
    const word& type,
    const dictionary& dict
)
{
    suspensionDepositionProperties_ = dict;
    coeffDict_ = suspensionDepositionProperties_.optionalSubDict(type + "Coeffs");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::suspensionDepositionModel::suspensionDepositionModel
(
    const word& type,
    const dictionary& suspensionDepositionProperties,
    const areaVectorField& Us,
    const areaScalarField& h,
    const areaScalarField& c,
    const areaVectorField& tau
)
:
    suspensionDepositionProperties_(suspensionDepositionProperties),
    coeffDict_
    (
        suspensionDepositionProperties_.optionalSubDict(type + "Coeffs")
    ),
    R_("R", dimless, suspensionDepositionProperties_),
    Us_(Us),
    h_(h),
    c_(c),
    tau_(tau),
    Sd_
    (
        IOobject
        (
            "Sd",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimVelocity)
    )
{
    Info<< "    with " << nl
        << "    " << R_ << endl;
}


// ************************************************************************* //
