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

Class
    Foam::suspensionDepositionModels::suspensionParkerFukushimaDeposition

Description
    Calculat the sedimentation rate in a turbidity current  following
    Parker et al., 1986

SourceFiles
    suspensionParkerFukushimaDeposition.C

Author
    Matthias Rauter matthias@rauter.it

\*---------------------------------------------------------------------------*/

#ifndef suspensionParkerFukushimaDeposition_H
#define suspensionParkerFukushimaDeposition_H

#include "suspensionDepositionModel.H"
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
namespace suspensionDepositionModels
{

/*---------------------------------------------------------------------------*\
                           Class suspensionParkerFukushimaDeposition Declaration
\*---------------------------------------------------------------------------*/

class suspensionParkerFukushimaDeposition
:
    public suspensionDepositionModel
{
    // Private data

        //- Form-factor for deposition curve
        dimensionedScalar Ds_;

        //- Form-factor for deposition curve
        dimensionedScalar nu_;

        //- Reference to the gravitation fields
        const areaVectorField& gs_;
        const areaScalarField& gn_;
        const areaScalarField& geff_;

public:

    //- Runtime type information
    TypeName("ParkerFukushimaDeposition");


    // Constructors

        //- Construct from components
        suspensionParkerFukushimaDeposition
        (
            const dictionary& frictionProperties,
            const areaVectorField& Us,
            const areaScalarField& h,
            const areaScalarField& c,
            const areaVectorField& tau
        );


    //- Destructor
    virtual ~suspensionParkerFukushimaDeposition() = default;


    // Member Functions

        //- Return Sink by deposition
        virtual const areaScalarField& Sd() const;

        //- Read depositionProperties dictionary
        virtual bool read(const dictionary& depositionProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frictionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
