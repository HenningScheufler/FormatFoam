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

Class
    Foam::functionObjects::totalVolume

Group
    grpUtilitiesFunctionObjects

Description
    calculate the total mass in a simulation. The mass is assumed to
    be the product of phase density c and flow height h. The phase
    density is optional.
Usage

SourceFiles
    totalVolume.C

\*---------------------------------------------------------------------------*/

#ifndef totalVolume_H
#define totalVolume_H

#include "regionFunctionObject.H"
#include "wordRes.H"
#include "Enum.H"
#include "gridfile.H"
#include "areaFieldsFwd.H"
#include "polyMesh.H"
#include "pointList.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class faMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                Class functionObjects::totalVolume Declaration
\*---------------------------------------------------------------------------*/

class totalVolume
:
    public regionFunctionObject
{
    // Private Data

        //- Reference to faMesh
        const faMesh& aMesh_;

        //- Name of the concentration field
        word cName_;

        //- Name of the flow height field
        word hName_;

        //- Output filename
        word tableOutput_;

        //- The file for results
        autoPtr<OFstream> dataFilePtr_;


    // Private Member Functions

        //- No copy construct
        totalVolume(const totalVolume&) = delete;

        //- No copy assignment
        void operator=(const totalVolume&) = delete;

public:

    //- Runtime type information
    TypeName("totalVolume");


    // Constructors

        //- Construct from Time and dictionary
        totalVolume
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~totalVolume() = default;


    // Member Functions

        //- Read the writeObjects data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write the registered objects
        virtual bool write();

        //- Update for changes of mesh - expires the surfaces
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for mesh point-motion - expires the surfaces
        virtual void movePoints(const polyMesh& mesh);

        //- Update for changes of mesh due to readUpdate - expires the surfaces
        virtual void readUpdate(const polyMesh::readUpdateState state);

        //- Mark intermediate surfaces and writers as needing an update.
        void expire();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
