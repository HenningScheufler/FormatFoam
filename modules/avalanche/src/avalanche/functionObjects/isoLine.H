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

Class
    Foam::functionObjects::isoLine

Group
    grpUtilitiesFunctionObjects

Description

Usage

SourceFiles
    isoLine.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_avalanche_isoLine_H
#define Foam_avalanche_isoLine_H

#include "regionFunctionObject.H"
#include "shapefile.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class faMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class functionObjects::isoLine Declaration
\*---------------------------------------------------------------------------*/

class isoLine
:
    public regionFunctionObject
{
    // Private Data

        //- Reference to faMesh
        const faMesh& aMesh_;

        //- outfile name
        word fileName_;

        //- field to be used for the isoline
        word fieldName_;

        //-iso-values to be used for the isoline
        List<scalar> isoValues_;

        //- Offset the outfile
        vector offset_;

public:

    //- No copy construct
    isoLine(const isoLine&) = delete;

    //- No copy assignment
    void operator=(const isoLine&) = delete;


    //- Runtime type information
    TypeName("isoLine");


    // Constructors

        //- Construct from Time and dictionary
        isoLine
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~isoLine() = default;


    // Member Functions

        //- Read the writeObjects data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write the registered objects
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
