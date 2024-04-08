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

Namespace
    Foam::suspensionDepositionModels

Description
    A namespace for various suspensionDeposition model implementations.

Class
    Foam::suspensionDepositionModel

Description
    An abstract base class for suspension deposition models

SourceFiles
    suspensionDepositionModel.C
    suspensionDepositionModelNew.C

Author
    Matthias Rauter matthias@rauter.it

\*---------------------------------------------------------------------------*/

#ifndef suspensionDepositionModel_H
#define suspensionDepositionModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"
#include "faMatrices.H"
#include "areaFieldsFwd.H"
#include "FieldFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class suspensionDepositionModel Declaration
\*---------------------------------------------------------------------------*/

class suspensionDepositionModel
{
protected:

    // Protected data

        dictionary suspensionDepositionProperties_;

        //- Model coefficients dictionary
        dictionary coeffDict_;

        //- The density ratio
        dimensionedScalar R_;

        //- Reference to the surface velocity field
        const areaVectorField& Us_;

        //- Reference to the flow height field
        const areaScalarField& h_;

        //- Reference to the sediment concentration field
        const areaScalarField& c_;

        //- Reference to the bottom stress field
        const areaVectorField& tau_;

        //- Sink term
        mutable areaScalarField Sd_;


    // Protected Member Functions

        //- Read/update the suspensionDepositionProperties and coeffDict dictionaries
        void readDict(const word& type, const dictionary& dict);


        //- Disallow copy construct
        suspensionDepositionModel(const suspensionDepositionModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const suspensionDepositionModel&) = delete;


public:

    //- Runtime type information
    TypeName("suspensionDepositionModel");


    // Declare run-time constructor selection table

#ifndef SWIG
        declareRunTimeSelectionTable
        (
            autoPtr,
            suspensionDepositionModel,
            dictionary,
            (
                const dictionary& suspensionDepositionProperties,
                const areaVectorField& Us,
                const areaScalarField& h,
                const areaScalarField& c,
                const areaVectorField& tau
            ),
            (suspensionDepositionProperties, Us, h, c, tau)
        );
#endif


    // Selectors

        //- Return a reference to the selected suspensionDeposition model
        static autoPtr<suspensionDepositionModel> New
        (
            const dictionary& suspensionDepositionProperties,
            const areaVectorField& Us,
            const areaScalarField& h,
            const areaScalarField& c,
            const areaVectorField& tau
        );


    // Constructors

        //- Construct from components
        suspensionDepositionModel
        (
            const word& type,
            const dictionary& suspensionDepositionProperties,
            const areaVectorField& Us,
            const areaScalarField& h,
            const areaScalarField& c,
            const areaVectorField& tau
        );


    //- Destructor
    virtual ~suspensionDepositionModel() = default;


    // Member Functions

        //- Read suspensionDepositionProperties dictionary
        virtual bool read(const dictionary& suspensionDepositionProperties) = 0;

        //- Return the friction properties dictionary
        const dictionary& suspensionDepositionProperties() const
        {
            return suspensionDepositionProperties_;
        }

        //- Const access to the model coefficients dictionary
        virtual const dictionary& coeffDict() const
        {
            return coeffDict_;
        }


        //- Return Sink by suspensionDeposition
        virtual const areaScalarField& Sd() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
