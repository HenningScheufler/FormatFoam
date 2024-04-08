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
    Foam::functionObjects::peakValues

Group
    grpUtilitiesFunctionObjects

Description
    Calculate the peak values of faFields over the runtime
Usage

SourceFiles
    peakValues.C

\*---------------------------------------------------------------------------*/

#ifndef avalance_peakValues_H
#define avalance_peakValues_H

#include "fvMeshFunctionObject.H"
#include "areaFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class faMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class functionObjects::peakValues Declaration
\*---------------------------------------------------------------------------*/

class peakValues
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Reference to faMesh
        const faMesh& aMesh_;

        //- The name given to the max field (output)
        word outputName_;

        //- Field for which peak values are calculated
        word fieldName_;


    // Private Member Functions

        //- No copy construct
        peakValues(const peakValues&) = delete;

        //- No copy assignment
        void operator=(const peakValues&) = delete;

public:

    //- Runtime type information
    TypeName("peakValues");


    // Constructors

        //- Construct from Time and dictionary
        peakValues
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~peakValues() = default;


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
