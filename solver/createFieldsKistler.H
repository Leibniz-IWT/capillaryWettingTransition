//The following two fields supports the Kistler contact angle model
volScalarField muEffKistler
(
    IOobject
    (
        "muEffKistler",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mixture.mu() + rho*turbulence->nut()
);

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Kistler contact angle model needs the sigma vol scalar field to be in the
// objectRegistry - even though sigma is already in the interface object
// below
volScalarField sigma
(
    IOobject
    (
        "sigmaKistler",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(transportProperties.lookup("sigma"))
);
