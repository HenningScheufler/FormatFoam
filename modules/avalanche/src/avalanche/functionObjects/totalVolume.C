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

\*---------------------------------------------------------------------------*/

#include "totalVolume.H"
#include "Time.H"
#include "polyMesh.H"
#include "faMesh.H"
#include "fvMesh.H"
#include "areaFields.H"
#include "ListOps.H"
#include "addToRunTimeSelectionTable.H"
#include "faCFD.H"
#include "HormannAgathos.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(totalVolume, 0);

    addRemovableToRunTimeSelectionTable
    (
        functionObject,
        totalVolume,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::totalVolume::totalVolume
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    #if (OPENFOAM < 2310)
    aMesh_(obr_.lookupObject<faMesh>("faMesh")),
    #else
    aMesh_(faMesh::mesh(refCast<const polyMesh>(obr_))),
    #endif
    cName_(),
    hName_("h"),
    tableOutput_("totalVolume.csv"),
    dataFilePtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::totalVolume::read(const dictionary& dict)
{
    regionFunctionObject::read(dict);

    cName_ = dict.getOrDefault<word>("cfield", word::null);
    hName_ = dict.getOrDefault<word>("hField", "h");

    tableOutput_ = dict.getOrDefault<word>("tableOutput", "totalVolume.csv");

    if (UPstream::master())
    {
        dataFilePtr_.reset(new OFstream(tableOutput_));
    }

    if (UPstream::master())
    {
        dataFilePtr_() << "time" << tab << "volume" << endl;
    }

    expire();

    return true;
}


bool Foam::functionObjects::totalVolume::execute()
{
    const auto& h =
        aMesh_.thisDb().lookupObject<areaScalarField>(hName_);

    const auto* concPtr =
        aMesh_.thisDb().cfindObject<areaScalarField>(cName_);

    scalar sumOfVolume = 0;

    if (concPtr)
    {
        const auto& c = *concPtr;
        sumOfVolume = gSum
        (
            h.primitiveField() * aMesh_.S()
          * c.primitiveField()
        );
    }
    else
    {
        sumOfVolume = gSum
        (
            h.primitiveField() * aMesh_.S()
        );
    }

    if (UPstream::master())
    {
        dataFilePtr_() << obr_.time().value()  << tab << sumOfVolume << endl;
    }

    return true;
}


bool Foam::functionObjects::totalVolume::write()
{
    return true;
}

void Foam::functionObjects::totalVolume::updateMesh(const mapPolyMesh& mpm)
{
    expire();
}

void Foam::functionObjects::totalVolume::movePoints(const polyMesh& mesh)
{
    expire();
}

void Foam::functionObjects::totalVolume::readUpdate(const polyMesh::readUpdateState state)
{
    expire();
}

void Foam::functionObjects::totalVolume::expire()
{

}


// ************************************************************************* //
