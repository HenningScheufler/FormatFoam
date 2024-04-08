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

Author
    Matthias Rauter matthias@rauter.it

\*---------------------------------------------------------------------------*/

#include "couplingModel.H"
#include "fvCFD.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(couplingModel, 0);
    defineRunTimeSelectionTable(couplingModel, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::couplingModel::readDict
(
    const word& type,
    const dictionary& dict
)
{
    couplingProperties_ = dict;
    coeffDict_ = couplingProperties_.optionalSubDict(type + "Coeffs");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::couplingModel::couplingModel
(
    const word& type,
    const dictionary& couplingProperties,
    const areaVectorField& Us1,
    const areaScalarField& h1,
    const areaScalarField& pb1,
    const areaVectorField& Us2,
    const areaScalarField& h2,
    const areaScalarField& c2
)
:
    couplingProperties_(couplingProperties),
    coeffDict_
    (
        couplingProperties_.optionalSubDict(type + "Coeffs")
    ),
    Us1_(Us1),
    h1_(h1),
    pb1_(pb1),
    Us2_(Us2),
    h2_(h2),
    c2_(c2),
    Sdp_
    (
        IOobject
        (
            "Sdp",
            Us1_.time().timeName(),
            Us1_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us1_.mesh(),
        dimensionedScalar(dimVelocity)
    )
{
    Info<< "    with " << endl;
}


// ************************************************************************* //
