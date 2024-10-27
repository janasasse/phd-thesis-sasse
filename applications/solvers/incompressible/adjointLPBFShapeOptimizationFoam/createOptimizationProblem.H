Info << "\nReading data for adjoint shape optimization\n" << endl;
IOdictionary optimizationProperties
(
    IOobject
    (
        "optimizationProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

word initialization(optimizationProperties.lookup("initialization"));
if (initialization == "yes")
{
    Info<< "Selecting run with shape optimisation turned off\n" << endl;
}
else
{
    Info<< "Selecting run with shape optimisation turned on\n" << endl;
}

Info<< "Reading initial geometry\n" << endl;
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading scaffold structure\n" << endl;
volScalarField alphaFixed
(
    IOobject
    (
        "alphaFixed",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

dimensionedScalar zeroSensitivity("0", dimVelocity*dimVelocity, 0.0);
dimensionedScalar alphaMax(optimizationProperties.lookup("alphaMax"));
dimensionedScalar zeroAlpha("0", dimless/dimTime, 0.0);

// sensitivity response step size * cell volume
dimensionedScalar lambdaU(optimizationProperties.lookup("lambdaU"));
dimensionedScalar lambdaT(optimizationProperties.lookup("lambdaT"));
dimensionedScalar lambdaC(optimizationProperties.lookup("lambdaC"));

// weights controlling optimisation objectives
dimensionedScalar wP(optimizationProperties.lookup("wP"));
dimensionedScalar wT(optimizationProperties.lookup("wT"));
dimensionedScalar wC(optimizationProperties.lookup("wC"));

// target values for optimisation objectives
dimensionedScalar targetTemperature(optimizationProperties.lookup("targetTemperature"));
dimensionedScalar targetConcentration(optimizationProperties.lookup("targetConcentration"));

if (wP.value() > 0.0)
{
    Info<< "Selecting optimisation objective 'minimal pressure loss'" << endl;
    Info<< "lambdaU: " << lambdaU.value()
        << "\n" << endl;
}

if (wT.value() > 0.0)
{
    Info<< "Selecting optimisation objective 'maximal thermal homogenisation'" << endl;
    Info<< "lambdaU: " << lambdaU.value()
        << ", lambdaT: " << lambdaT.value()
        << ", target temperature: " << targetTemperature.value()
        << " K\n" << endl;
}

if (wC.value() > 0.0)
{
    Info<< "Selecting optimisation objective 'maximal material homogenisation'" << endl;
    Info<< "lambdaU: " << lambdaU.value()
        << ", lambdaC: " << lambdaC.value()
        << ", target concentration: " << targetConcentration.value()
        << " m^-3\n" << endl;
}

Info<< "Creating Sensitivitiy Field\n" << endl;
volScalarField sens
(
    IOobject
    (
        "sensitivity",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    1e-5 * (lambdaU*(Ua&U))
);

Info<< "Creating fields for geometric constraints\n" << endl;

volScalarField alphaDiff
(
    IOobject
    (
        "alphaDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alphaDiff", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

volScalarField overhangCriterion
(
    IOobject
    (
        "overhangCriterion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("overhangCriterion", dimensionSet(0,0,0,0,0,0,0), 0.0)
);

volScalarField minLengthCriterion
(
    IOobject
    (
        "minLengthCriterion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("minLengthCriterion", dimensionSet(0,0,0,0,0,0,0), 0.0)
);

volScalarField rho_physical
(
    IOobject
    (
        "rho_physical",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rho_physical", dimensionSet(0,0,0,0,0,0,0), 0.0)
);

volScalarField rho_im
(
    IOobject
    (
        "rho_im",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rho_im", dimensionSet(0,0,0,0,0,0,0), 0.0)
);
