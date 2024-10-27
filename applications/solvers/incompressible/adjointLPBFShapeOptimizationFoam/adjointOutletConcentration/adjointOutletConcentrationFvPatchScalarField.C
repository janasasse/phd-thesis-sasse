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

#include "adjointOutletConcentrationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointOutletConcentrationFvPatchScalarField::
adjointOutletConcentrationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::adjointOutletConcentrationFvPatchScalarField::
adjointOutletConcentrationFvPatchScalarField
(
    const adjointOutletConcentrationFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::adjointOutletConcentrationFvPatchScalarField::
adjointOutletConcentrationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


Foam::adjointOutletConcentrationFvPatchScalarField::
adjointOutletConcentrationFvPatchScalarField
(
    const adjointOutletConcentrationFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointOutletConcentrationFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U");

    const fvPatchField<scalar>& Concp =
        patch().lookupPatchField<volScalarField, scalar>("Conc");

    const fvPatchField<scalar>& Concap =
        patch().lookupPatchField<volScalarField, scalar>("Conca");

    scalarField Un = Up & patch().nf();

    const dictionary& optimizationProperties = db().lookupObject<IOdictionary>("optimizationProperties");
    dimensionedScalar wC(optimizationProperties.lookup("wC"));
    dimensionedScalar targetConcentration(optimizationProperties.lookup("targetConcentration"));

    scalarField costFunctionConcentration = (Concp - targetConcentration.value());

    // neglect second order term to improve stability
    scalarField tmp=((Un * Concap) + wC.value() * costFunctionConcentration);

    operator==(tmp);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::adjointOutletConcentrationFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        adjointOutletConcentrationFvPatchScalarField
    );
}

// ************************************************************************* //
