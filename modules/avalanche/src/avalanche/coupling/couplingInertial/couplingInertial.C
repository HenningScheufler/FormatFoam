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

#include "fvCFD.H"
#include "faCFD.H"
#include "couplingInertial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace couplingModels
{
    defineTypeNameAndDebug(couplingInertial, 0);
    addToRunTimeSelectionTable(couplingModel, couplingInertial, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::couplingModels::couplingInertial::couplingInertial
(
    const dictionary& couplingProperties,
    const areaVectorField& Us1,
    const areaScalarField& h1,
    const areaScalarField& pb1,
    const areaVectorField& Us2,
    const areaScalarField& h2,
    const areaScalarField& c2
)
:
    couplingModel(type(), couplingProperties, Us1, h1, pb1, Us2, h2, c2),
    I0_("I0", dimless, coeffDict_),
    u0_("u0", dimless, coeffDict_),
    d_("d", dimLength, coeffDict_),
    rhos_("rhos", dimDensity, coeffDict_),
    I_
    (
        IOobject
        (
            "I",
            Us1_.time().timeName(),
            Us1_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Us1_.mesh(),
        dimensionedScalar(dimless)
    )
{
    Info << "    " << u0_ << nl 
         << "    " << d_ << nl 
         << "    " << rhos_ << nl << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::areaScalarField&
Foam::couplingModels::couplingInertial::Sdp() const
{
    //Inertial number at the base, assuming it is constant along flow depth
    I_ = 5./2.*mag(Us1_)/(h1_+dimensionedScalar(dimLength, 1e-2))*d_*sqrt(rhos_/(pb1_+dimensionedScalar(dimPressure, 1e-2)));
    Sdp_ = mag(Us1_)*u0_*mag(I_-I0_)*h1_/(h1_+dimensionedScalar(dimLength, 1e-2));
    return Sdp_;
}


bool Foam::couplingModels::couplingInertial::read
(
    const dictionary& couplingProperties
)
{
    readDict(type(), couplingProperties);

    coeffDict_.readEntry("u0", u0_);
    coeffDict_.readEntry("I0", I0_);
    coeffDict_.readEntry("d", d_);
    coeffDict_.readEntry("rhos", rhos_);

    return true;
}


// ************************************************************************* //
