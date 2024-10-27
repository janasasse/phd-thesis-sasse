Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading temperature field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading concentration field Conc\n" << endl;
volScalarField Conc
(
    IOobject
    (
        "Conc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volTensorField tau
(
    IOobject
    (
        "tau",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("tau", dimensionSet(0,2,-2,0,0,0,0), tensor::zero)
);

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Reading data for fluid and solid specific heat capacity c\n" <<endl;
    dimensionedScalar cfluid
    (
    transportProperties.lookup("cfluid")
    );

    dimensionedScalar csolid
    (
    transportProperties.lookup("csolid")
    );

    volScalarField c
    (
        IOobject
        (
            "c",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("c", dimensionSet(0,2,-2,-1,0,0,0), cfluid.value())
    );

Info<< "Reading data for fluid and solid thermal diffusivities DT\n" <<endl;
    dimensionedScalar DTfluid
    (
    transportProperties.lookup("DTfluid")
    );

    dimensionedScalar DTsolid
    (
    transportProperties.lookup("DTsolid")
    );

    volScalarField DT
    (
        IOobject
        (
            "DTfield",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("DTfield", dimensionSet(0,2,-1,0,0,0,0), DTfluid.value())
    );

Info<< "Reading data for fluid and solid diffusion coefficients DConc\n" <<endl;
    dimensionedScalar DConcfluid
    (
    transportProperties.lookup("DConcfluid")
    );

    dimensionedScalar DConcsolid
    (
    transportProperties.lookup("DConcsolid")
    );

    dimensionedScalar DConca
    (
    transportProperties.lookup("DConca")
    );

    volScalarField DConc
    (
        IOobject
        (
            "DConcfield",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("DConcfield", dimensionSet(0,2,-1,0,0,0,0), DConcfluid.value())
    );

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


Info<< "Reading field pa\n" << endl;
volScalarField pa
(
    IOobject
    (
        "pa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Ua\n" << endl;
volVectorField Ua
(
    IOobject
    (
        "Ua",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading temperature field Ta\n" << endl;
volScalarField Ta
(
    IOobject
    (
        "Ta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading concentration field Conca\n" << endl;
volScalarField Conca
(
    IOobject
    (
        "Conca",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhia.H"


label paRefCell = 0;
scalar paRefValue = 0.0;
setRefCell
(
    pa,
    simple.dict(),
    paRefCell,
    paRefValue
);
mesh.setFluxRequired(pa.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createFvOptions.H"

#include "createOptimizationProblem.H"

const labelList& inletCells = mesh.boundary()["inlet"].faceCells();
//const labelList& outletCells = mesh.boundary()["outlet"].faceCells();

zeroCells(alpha, inletCells);
//zeroCells(alpha, outletCells);


