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

\*---------------------------------------------------------------------------*/

#include "adjointOutletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}


Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const adjointOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const adjointOutletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::adjointOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U");

    const fvPatchField<vector>& Uap =
        patch().lookupPatchField<volVectorField, vector>("Ua");

    const fvPatchField<scalar>& Tp =
        patch().lookupPatchField<volScalarField, scalar>("T");

    const fvPatchField<scalar>& Tap =
        patch().lookupPatchField<volScalarField, scalar>("Ta");

    const fvPatchField<scalar>& Concp =
        patch().lookupPatchField<volScalarField, scalar>("Conc");

    const fvPatchField<scalar>& Concap =
        patch().lookupPatchField<volScalarField, scalar>("Conca");

    vectorField Ut = Up - patch().nf() * (Up&patch().nf());
    tmp<vectorField> tnf(patch().nf());
    const vectorField& nf = tnf();
    tmp<vectorField> tUan(nf * (Uap&nf));
    const vectorField& Uan = tUan();
    vectorField Uat = (Uap - Uan);
    scalarField Uans = (Uap & patch().nf());

    vectorField UanUt = (Uans*Ut);
    vectorField UanUat = (Uans*Uat);

    // neglect second order term to improve stability
    vectorField::operator == (UanUt - UanUat + (Tp * Tap * patch().nf()) + (Concp * Concap * patch().nf()));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::adjointOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointOutletVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
