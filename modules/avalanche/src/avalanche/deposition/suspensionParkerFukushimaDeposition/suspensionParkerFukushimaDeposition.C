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

#include "fvCFD.H"
#include "faCFD.H"
#include "suspensionParkerFukushimaDeposition.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace suspensionDepositionModels
{
    defineTypeNameAndDebug(suspensionParkerFukushimaDeposition, 0);
    addToRunTimeSelectionTable(suspensionDepositionModel, suspensionParkerFukushimaDeposition, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::suspensionDepositionModels::suspensionParkerFukushimaDeposition::suspensionParkerFukushimaDeposition
(
    const dictionary& depositionProperties,
    const areaVectorField& Us,
    const areaScalarField& h,
    const areaScalarField& c,
    const areaVectorField& tau
)
:
    suspensionDepositionModel(type(), depositionProperties, Us, h, c, tau),
    Ds_("Ds", coeffDict_),
    nu_("nu", coeffDict_),
    gs_(Us.db().lookupObject<areaVectorField>("gs")),
    gn_(Us.db().lookupObject<areaScalarField>("gn")),
    geff_(Us.db().lookupObject<areaScalarField>("geff"))
{
    Info<< "    " << Ds_ << nl
        << "    " << nu_ << nl
        << "    " << Ds_ << nl << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::areaScalarField&
Foam::suspensionDepositionModels::suspensionParkerFukushimaDeposition::Sd() const
{

    const areaScalarField vs(R_*geff_*Ds_*Ds_/18./nu_); //Terminal velocity
    const areaScalarField mu(sqrt(mag(tau_))/(vs+dimensionedScalar(dimVelocity, SMALL))); //Parker et al. (1986), Eq. (21)
    const areaScalarField r0(1.+31.5*pow(mu+SMALL, -1.46)); //Parker et al. (1986), Eq. (20)
    //See Parker et al. (1986), Eq. (4), last term
    Sd_ = vs*r0*c_;

    return Sd_;
}


bool Foam::suspensionDepositionModels::suspensionParkerFukushimaDeposition::read
(
    const dictionary& depositionProperties
)
{
    readDict(type(), depositionProperties);

    coeffDict_.readEntry("R", R_);
    coeffDict_.readEntry("Ds", Ds_);
    coeffDict_.readEntry("nu", nu_);

    return true;
}


// ************************************************************************* //
