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

\*---------------------------------------------------------------------------*/

#include "peakValues.H"
#include "Time.H"
#include "areaFields.H"
#include "ListOps.H"
#include "addToRunTimeSelectionTable.H"
#include "shapefile.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(peakValues, 0);

    addRemovableToRunTimeSelectionTable
    (
        functionObject,
        peakValues,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::peakValues::peakValues
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    #if (OPENFOAM < 2310)
    aMesh_(mesh_.lookupObject<faMesh>("faMesh")),
    #else
    aMesh_(faMesh::mesh(mesh_)),
    #endif
    outputName_("max_"),
    fieldName_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::peakValues::read(const dictionary& dict)
{
    regionFunctionObject::read(dict);

    dict.readEntry("field", fieldName_);
    dict.readEntry("output", outputName_);

    return true;
}


bool Foam::functionObjects::peakValues::execute()
{
    // scalar
    {
        const auto* phiPtr =
            aMesh_.thisDb().cfindObject<areaScalarField>(fieldName_);

        if (phiPtr)
        {
            const areaScalarField& phi = *phiPtr;

            auto* maxField =
                aMesh_.thisDb().getObjectPtr<areaScalarField>(outputName_);

            if (!maxField)
            {
                maxField = new areaScalarField
                (
                    IOobject
                    (
                        outputName_,
                        obr_.time().timeName(),
                        aMesh_.thisDb(),
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE,
                        IOobject::REGISTER
                    ),
                    phi
                );
                regIOobject::store(maxField);

                Info<< "Initializing peak Value field \"" << outputName_
                    << "\" for areaScalarField \"" << fieldName_ << "\"" << nl;
            }
            else
            {
                *maxField = max(phi, *maxField);
            }

            return true;
        }
    }

    // vector
    {
        const auto* phiPtr =
            aMesh_.thisDb().cfindObject<areaVectorField>(fieldName_);

        if (phiPtr)
        {
            const auto& phi = *phiPtr;

            auto* maxField =
                aMesh_.thisDb().getObjectPtr<areaScalarField>(outputName_);

            if (!maxField)
            {
                maxField = new areaScalarField
                (
                    IOobject
                    (
                        outputName_,
                        obr_.time().timeName(),
                        aMesh_.thisDb(),
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE,
                        IOobject::REGISTER
                    ),
                    mag(phi)
                );

                regIOobject::store(maxField);

                Info << "Initializing peak Value field \"" << outputName_
                    << "\" for areaVectorField \"" << fieldName_ << "\"" << nl;
            }
            else
            {
                *maxField = max(mag(phi), *maxField);
            }
        }

        return true;
    }

    // Fallthrough
    return false;
}


bool Foam::functionObjects::peakValues::write()
{
    return true;
}


// ************************************************************************* //
