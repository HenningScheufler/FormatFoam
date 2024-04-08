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
    Foam::couplingModels::couplingInertial
Description
    A coupling model for mixed snow avalanches following the scaling with the
    Inertial number I.

    Publication in preparation.

SourceFiles
    couplingInertial.C

Author
    Matthias Rauter matthias@rauter.it

\*---------------------------------------------------------------------------*/

#ifndef couplingInertial_H
#define couplingInertial_H

#include "couplingModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"
#include "faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace couplingModels
{

/*---------------------------------------------------------------------------*\
                           Class Erosionenergy Declaration
\*---------------------------------------------------------------------------*/

class couplingInertial
:
    public couplingModel
{
    // Private data


        //- reference inertial number
        dimensionedScalar I0_;

        //- coupling volume flux at I-I0 = 1
        dimensionedScalar u0_;

        //- particle diameter for calculation I
        dimensionedScalar d_;

        //- rhos for calculating I
        dimensionedScalar rhos_;

        //- Inertial number
        mutable areaScalarField I_;
public:

    //- Runtime type information
    TypeName("couplingInertial");


    // Constructors

        //- Construct from components
        couplingInertial
        (
            const dictionary& couplingProperties,
            const areaVectorField& Us1,
            const areaScalarField& h1,
            const areaScalarField& pb1,
            const areaVectorField& Us2,
            const areaScalarField& h2,
            const areaScalarField& c2
        );


    //- Destructor
    virtual ~couplingInertial() = default;


    // Member Functions

        //- Return the volume flux
        virtual const areaScalarField& Sdp() const;

        //- Read couplingProperties dictionary
        virtual bool read(const dictionary& couplingProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frictionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
