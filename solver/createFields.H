#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
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

#include "createPhi.H"


//The following two fields supports the Kistler contact angle model
volScalarField muEffKistler
(
    IOobject
    (
        "muEffKistler",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
 //       IOobject::NO_WRITE
        IOobject::AUTO_WRITE
    ),
    mesh
    //,
    //mixture.mu() + rho*turbulence->nut()
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

    

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);
volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());
const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();



// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

    // Construct interface from alpha1 distribution
    interfaceSigmaProperties interfaceSigma(alpha1, U, transportProperties);

    //This dictionary contains important props for phase change:
    const dictionary& surfaceTensionForceProperties
    (
        transportProperties.subDict("surfaceTensionForce")
    );


    //Initialize the surface tension force model:
    autoPtr<surfaceTensionForceModel> stfModel
    (
        surfaceTensionForceModel::New
        (
            surfaceTensionForceProperties.lookup("model"),
            surfaceTensionForceProperties,
            interfaceSigma,
            alpha1
        )
    );

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
//#include "createFieldsKistler.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

//Creating field for gradU
    Info<< "Reading field gradU\n" << endl;
	volTensorField gradU
	(
	    IOobject
	    (
		"gradU",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,	//is read, needs to be in folder 0
		IOobject::AUTO_WRITE	//is written
	    ),
	    mesh
	);

    //Creating field for D
    Info<< "Reading field D\n" << endl;
	volTensorField D
	(
	    IOobject
	    (
		"D",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,	//is read, needs to be in folder 0
		IOobject::AUTO_WRITE	//is written
	    ),
	    mesh
	);

    //Creating field for shearRate
    Info<< "Reading field shearRate\n" << endl;
	volVectorField shearRate
	(
	    IOobject
	    (
		"shearRate",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,	//is read, needs to be in folder 0
		IOobject::AUTO_WRITE	//is written
	    ),
	    mesh
	);

    //Creating field for tauShear
    Info<< "Reading field tauShear\n" << endl;
	volVectorField tauShear
	(
	    IOobject
	    (
		"tauShear",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,	//is read, needs to be in folder 0
		IOobject::AUTO_WRITE	//is written
	    ),
	    mesh
	);

    //Creating field for strainRate
    Info<< "Reading field strainRate\n" << endl;
	volVectorField strainRate
	(
	    IOobject
	    (
		"strainRate",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,	//is read, needs to be in folder 0
		IOobject::AUTO_WRITE	//is written
	    ),
	    mesh
	);

    //Creating field for tauStrain
    Info<< "Reading field tauStrain\n" << endl;
	volVectorField tauStrain
	(
	    IOobject
	    (
		"tauStrain",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,	//is read, needs to be in folder 0
		IOobject::AUTO_WRITE	//is written
	    ),
	    mesh
	);


#include "createMRF.H"
#include "createFvOptions.H"
#include "createIsoAdvection.H"
