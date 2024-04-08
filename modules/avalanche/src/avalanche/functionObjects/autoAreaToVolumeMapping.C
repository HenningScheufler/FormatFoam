/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | avalanche module
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Matthias Rauter
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "autoAreaToVolumeMapping.H"
#include "Time.H"
#include "faCFD.H"
#include "fvCFD.H"
#include "polyMesh.H"
#include "faMesh.H"
#include "areaFields.H"
#include "ListOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(autoAreaToVolumeMapping, 0);

    addRemovableToRunTimeSelectionTable
    (
        functionObject,
        autoAreaToVolumeMapping,
        dictionary
    );
}
}

const Foam::Enum
<
    Foam::functionObjects::autoAreaToVolumeMapping::writeOption
>
Foam::functionObjects::autoAreaToVolumeMapping::writeOptionNames_
({
    { writeOption::AUTO_WRITE, "autoWrite" },
    { writeOption::NO_WRITE, "noWrite" },
    { writeOption::ANY_WRITE, "anyWrite" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::autoAreaToVolumeMapping::autoAreaToVolumeMapping
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
    writeOption_(ANY_WRITE),
    objectNames_(),
    prefix_("fa_"),
    vsm_(aMesh_)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::autoAreaToVolumeMapping::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (dict.found("field"))
    {
        objectNames_.resize(1);
        dict.readEntry("field", objectNames_.first());
    }
    else if (dict.found("fields"))
    {
        dict.readEntry("fields", objectNames_);
    }
    else
    {
        dict.readEntry("objects", objectNames_);
    }

    writeOption_ = writeOptionNames_.getOrDefault
    (
        "writeOption",
        dict,
        writeOption::ANY_WRITE
    );

    dict.readIfPresent("prefix", prefix_);

    return true;
}


bool Foam::functionObjects::autoAreaToVolumeMapping::execute()
{
    return true;
}


bool Foam::functionObjects::autoAreaToVolumeMapping::write()
{
    Info << type() << " " << name() << " write:" << nl;

    if (!obr_.time().writeTime())
    {
        return false;
    }

    // Get selection
    const wordList selectedNames
    (
        aMesh_.thisDb().sortedNames(objectNames_)
    );

    // Warning if anything was missed
    bitSet missed(objectNames_.size());

    label index = 0;
    for (const wordRe& select : objectNames_)
    {
        if (!ListOps::found(selectedNames, select))
        {
            missed.set(index);
        }
        ++index;
    }

    if (missed.any())
    {
        WarningInFunction
            << "No corresponding selection for "
            << flatOutput(subset(missed, objectNames_)) << nl
            << "Available objects in database:" << nl
            << aMesh_.thisDb().sortedToc()
            << endl;
    }

    for (const word& objName : selectedNames)
    {
        const regIOobject& obj =
            aMesh_.thisDb().lookupObject<regIOobject>(objName);

        if (!fieldTypes::area.contains(obj.type()))
        {
            continue;
        }

        switch (writeOption_)
        {
            case AUTO_WRITE:
            {
                if (obj.writeOpt() != IOobject::AUTO_WRITE)
                {
                    continue;
                }
                break;
            }
            case NO_WRITE:
            {
                if (obj.writeOpt() != IOobject::NO_WRITE)
                {
                    continue;
                }
                break;
            }
            case ANY_WRITE:
            {
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown writeOption "
                    << writeOptionNames_[writeOption_]
                    << ". Valid writeOption types are "
                    << writeOptionNames_
                    << exit(FatalError);

                continue;
                break;
            }
        }
        if (obj.name().ends_with("_0"))
        {
            Info<< "    ignoring old " << obj.type()
                << ' ' << obj.name() << endl;
            continue;
        }
        else if
        (
            obj.writeOpt() == IOobject::AUTO_WRITE
         && obr_.time().writeTime()
        )
        {
            Info<< "    automatically writing " << obj.type()
                << ' ' << obj.name() << endl;
        }
        else
        {
            Info<< "    writing " << obj.type()
                << ' ' << obj.name() << endl;
        }

        #undef doLocalCode
        #define doLocalCode(Type)                                             \
        {                                                                     \
            const auto* afld = isA<AreaField<Type>>(obj);                     \
                                                                              \
            if (afld)                                                         \
            {                                                                 \
                VolumeField<Type> vfld                                        \
                (                                                             \
                    IOobject                                                  \
                    (                                                         \
                        prefix_ + obj.name(),                                 \
                        mesh_.time().timeName(),                              \
                        mesh_.thisDb(),                                       \
                        IOobject::NO_READ,                                    \
                        IOobject::NO_WRITE,                                   \
                        IOobject::NO_REGISTER                                 \
                    ),                                                        \
                    mesh_,                                                    \
                    Zero,                                                     \
                    afld->dimensions()                                        \
                );                                                            \
                                                                              \
                vsm_.mapToVolume(*afld, vfld.boundaryFieldRef());             \
                vfld.write();                                                 \
                continue;                                                     \
            }                                                                 \
        }

        doLocalCode(scalar);
        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);

        #undef doLocalCode
    }

    Info << endl;

    return true;
}


// ************************************************************************* //
