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
#include "ambientParkerFukushimaEntrainment.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ambientEntrainmentModels
{
    defineTypeNameAndDebug(ambientParkerFukushimaEntrainment, 0);
    addToRunTimeSelectionTable(ambientEntrainmentModel, ambientParkerFukushimaEntrainment, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ambientEntrainmentModels::ambientParkerFukushimaEntrainment::ambientParkerFukushimaEntrainment
(
    const dictionary& entrainmentProperties,
    const areaVectorField& Us,
    const areaScalarField& h,
    const areaScalarField& c
)
:
    ambientEntrainmentModel(type(), entrainmentProperties, Us, h, c),
    ewf_("ewf", dimless, coeffDict_),
    Ri0_("Ri0", dimless, coeffDict_),
    geff_(Us_.db().lookupObject<areaScalarField>("geff"))
{
    Info << "    " << ewf_ << nl
         << "    " << Ri0_ << nl << endl;

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::areaScalarField&
Foam::ambientEntrainmentModels::ambientParkerFukushimaEntrainment::Sm() const
{
    dimensionedScalar u0(dimVelocity, 1e-5);
    const areaScalarField Ri(R_*geff_*c_*h_/(magSqr(Us_)+sqr(u0))); //Parker et al. (1986) Eq. (10)
    Sm_ = ewf_/(Ri0_+Ri)*mag(Us_); //Parker et al. (1986) Eq. (19)

    return Sm_;
}


bool Foam::ambientEntrainmentModels::ambientParkerFukushimaEntrainment::read
(
    const dictionary& entrainmentProperties
)
{
    readDict(type(), entrainmentProperties);

    coeffDict_.readEntry("R", R_);
    coeffDict_.readEntry("ewf", ewf_);
    coeffDict_.readEntry("Ri0", Ri0_);

    return true;
}


// ************************************************************************* //
