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

Author
    Matthias Rauter matthias@rauter.it

\*---------------------------------------------------------------------------*/

#include "ambientEntrainmentModel.H"
#include "fvCFD.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ambientEntrainmentModel, 0);
    defineRunTimeSelectionTable(ambientEntrainmentModel, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::ambientEntrainmentModel::readDict
(
    const word& type,
    const dictionary& dict
)
{
    ambientEntrainmentProperties_ = dict;
    coeffDict_ = ambientEntrainmentProperties_.optionalSubDict(type + "Coeffs");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ambientEntrainmentModel::ambientEntrainmentModel
(
    const word& type,
    const dictionary& ambientEntrainmentProperties,
    const areaVectorField& Us,
    const areaScalarField& h,
    const areaScalarField& c
)
:
    ambientEntrainmentProperties_(ambientEntrainmentProperties),
    coeffDict_
    (
        ambientEntrainmentProperties_.optionalSubDict(type + "Coeffs")
    ),
    R_("R", dimless, ambientEntrainmentProperties_),
    Us_(Us),
    h_(h),
    c_(c),
    Sm_
    (
        IOobject
        (
            "Sm",
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
        << "    " << R_  << nl;
}


// ************************************************************************* //
