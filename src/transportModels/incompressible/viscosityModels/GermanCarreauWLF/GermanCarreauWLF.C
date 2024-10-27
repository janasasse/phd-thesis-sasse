/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "GermanCarreauWLF.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(GermanCarreauWLF, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        GermanCarreauWLF,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::GermanCarreauWLF::calcNu() const
{
    const volScalarField& Tcurr = U_.mesh().lookupObject<volScalarField>("T");

    const dimensionedScalar C1("C1", dimensionSet(0,0,0,0,0,0,0), 8.86);
    const dimensionedScalar C2("C2", dimensionSet(0,0,0,1,0,0,0), 101.6);

    volScalarField tmpTs
    (
        IOobject
        (
            "tmpTs",
            U_.time().timeName(),
			U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        Ts_
    );

    volScalarField tmpTm
    (
        IOobject
        (
            "tmpTm",
            U_.time().timeName(),
			U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        Tm_
    );

    volScalarField termB = (C1*(Tcurr - tmpTs)) / (C2+(Tcurr - tmpTs));
    volScalarField termA = (C1*(tmpTm - tmpTs)) / (C2+(tmpTm - tmpTs));
    volScalarField alpha_shift = pow(scalar(10), (termA - termB));

    return 
        (nuA_*alpha_shift) 
      / pow(scalar(1) + alpha_shift*B_*strainRate(), C_);
	
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::GermanCarreauWLF::GermanCarreauWLF
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    GermanCarreauWLFCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    nuA_("nuA", dimViscosity, GermanCarreauWLFCoeffs_),
    B_("B", dimTime, GermanCarreauWLFCoeffs_),
    C_("C", dimless, GermanCarreauWLFCoeffs_),
    Ts_("Ts", dimTemperature, GermanCarreauWLFCoeffs_),
    Tm_("Tm", dimTemperature, GermanCarreauWLFCoeffs_),

    nu_
	(
        IOobject
		(
			name,
			U_.time().timeName(),
			U_.db(),
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		
		),	
       	calcNu()
	)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::GermanCarreauWLF::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    GermanCarreauWLFCoeffs_.lookup("nuA")   >> nuA_;
    GermanCarreauWLFCoeffs_.lookup("B")     >> B_;
    GermanCarreauWLFCoeffs_.lookup("C")     >> C_;
    GermanCarreauWLFCoeffs_.lookup("Ts")    >> Ts_;
    GermanCarreauWLFCoeffs_.lookup("Tm")    >> Tm_;

    return true;
}


// ************************************************************************* //
