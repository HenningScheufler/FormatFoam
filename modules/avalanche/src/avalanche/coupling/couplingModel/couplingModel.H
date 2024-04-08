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
    Foam::couplingModel

Description
    A namespace for various coupling model implementations for mixed snow
    avalanches. This model couples two shallow flow layers.

Class
    Foam::couplingModel

Description
    An abstract base class for coupling models.

SourceFiles
    couplingModel.C
    couplingModelNew.C

Author
    Matthias Rauter matthias@rauter.it

\*---------------------------------------------------------------------------*/

#ifndef couplingModel_H
#define couplingModel_H

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
                         Class couplingModel Declaration
\*---------------------------------------------------------------------------*/

class couplingModel
{
protected:

    // Protected data

        dictionary couplingProperties_;

        //- Model coefficients dictionary
        dictionary coeffDict_;

        //- Reference to the surface velocity field
        const areaVectorField& Us1_;

        //- Reference to the flow height field
        const areaScalarField& h1_;

        //- Reference to the basal pressure field
        const areaScalarField& pb1_;

        //- Reference to the surface velocity field
        const areaVectorField& Us2_;

        //- Reference to the flow height field
        const areaScalarField& h2_;

        //- Reference to the phase fraction field
        const areaScalarField& c2_;

        //- Flux term
        mutable areaScalarField Sdp_;


    // Protected Member Functions

        //- Read/update the couplingProperties and coeffDict dictionaries
        void readDict(const word& type, const dictionary& dict);


        //- Disallow copy construct
        couplingModel(const couplingModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const couplingModel&) = delete;


public:

    //- Runtime type information
    TypeName("couplingModel");


    // Declare run-time constructor selection table

#ifndef SWIG
        declareRunTimeSelectionTable
        (
            autoPtr,
            couplingModel,
            dictionary,
            (
                const dictionary& couplingProperties,
                const areaVectorField& Us1,
                const areaScalarField& h1,
                const areaScalarField& pb1,
                const areaVectorField& Us2,
                const areaScalarField& h2,
                const areaScalarField& c2
            ),
            (couplingProperties, Us1, h1, pb1, Us2, h2, c2)
        );
#endif


    // Selectors

        //- Return a reference to the selected coupling model
        static autoPtr<couplingModel> New
        (
            const dictionary& couplingProperties,
            const areaVectorField& Us1,
            const areaScalarField& h1,
            const areaScalarField& pb1,
            const areaVectorField& Us2,
            const areaScalarField& h2,
            const areaScalarField& c2
        );


    // Constructors

        //- Construct from components
        couplingModel
        (
            const word& type,
            const dictionary& couplingProperties,
            const areaVectorField& Us1,
            const areaScalarField& h1,
            const areaScalarField& pb1,
            const areaVectorField& Us2,
            const areaScalarField& h2,
            const areaScalarField& c2
        );


    //- Destructor
    virtual ~couplingModel() = default;


    // Member Functions

        //- Read couplingProperties dictionary
        virtual bool read(const dictionary& couplingProperties) = 0;

        //- Return the coupling properties dictionary
        const dictionary& couplingProperties() const
        {
            return couplingProperties_;
        }

        //- Const access to the model coefficients dictionary
        virtual const dictionary& coeffDict() const
        {
            return coeffDict_;
        }

        //- Return the Source by coupling
        virtual const areaScalarField& Sdp() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
