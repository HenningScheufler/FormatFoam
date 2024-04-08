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
#include "suspensionEntrainmentOff.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace suspensionEntrainmentModels
{
    defineTypeNameAndDebug(suspensionEntrainmentOff, 0);
    addToRunTimeSelectionTable(suspensionEntrainmentModel, suspensionEntrainmentOff, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::suspensionEntrainmentModels::suspensionEntrainmentOff::suspensionEntrainmentOff
(
    const dictionary& entrainmentProperties,
    const areaVectorField& Us,
    const areaScalarField& h,
    const areaScalarField& hentrain,
    const areaScalarField& c,
    const areaVectorField& tau
)
:
    suspensionEntrainmentModel(type(), entrainmentProperties, Us, h, hentrain, c, tau)
{
    Info << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::areaScalarField&
Foam::suspensionEntrainmentModels::suspensionEntrainmentOff::Sm() const
{
    return Sm_;
}


bool Foam::suspensionEntrainmentModels::suspensionEntrainmentOff::read
(
    const dictionary& entrainmentProperties
)
{
    readDict(type(), entrainmentProperties);

    return true;
}


// ************************************************************************* //
