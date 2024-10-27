/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::adjointOutletConcentrationFvPatchScalarField

Description

SourceFiles
    adjointOutletConcentrationFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef adjointOutletConcentrationFvPatchScalarField_H
#define adjointOutletConcentrationFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class adjointOutletConcentrationFvPatch Declaration
\*---------------------------------------------------------------------------*/

class adjointOutletConcentrationFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("adjointOutletConcentration");


    // Constructors

        //- Construct from patch and internal field
        adjointOutletConcentrationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        adjointOutletConcentrationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given adjointOutletConcentrationFvPatchScalarField
        //  onto a new patch
        adjointOutletConcentrationFvPatchScalarField
        (
            const adjointOutletConcentrationFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new adjointOutletConcentrationFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        adjointOutletConcentrationFvPatchScalarField
        (
            const adjointOutletConcentrationFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new adjointOutletConcentrationFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
