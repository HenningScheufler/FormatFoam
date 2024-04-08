    Info << "Reading field h2" << endl;
    areaScalarField h2
    (
        IOobject
        (
            "h2",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh
    );


    Info << "Reading field Us2" << endl;
    areaVectorField Us2
    (
        IOobject
        (
            "Us2",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh
    );

    areaVectorField UsDiff
    (
        IOobject
        (
            "UsDiff",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us2-Us1
    );

    Info << "Reading field c2" << endl;
    areaScalarField c2
    (
        IOobject
        (
            "c2",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh
    );

    areaScalarField powderLayerAlone
    (
        IOobject
        (
            "powderLayerAlone",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh,
        dimensionedScalar(dimless)
    );

    edgeScalarField phis2
    (
        IOobject
        (
            "phis2",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fac::interpolate(Us2) & aMesh.Le()
    );

    edgeScalarField phi2s2
    (
        IOobject
        (
            "phi2s2",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fac::interpolate(h2*Us2) & aMesh.Le()
    );

    areaScalarField geff
    (
        IOobject
        (
            "geff",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        gn
    );

    areaScalarField Sdp
    (
        IOobject
        (
            "Sdp",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh,
        dimensionedScalar(dimVelocity)
    );

    areaScalarField Spd
    (
        IOobject
        (
            "Spd",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh,
        dimensionedScalar(dimVelocity)
    );

    areaScalarField Sap
    (
        IOobject
        (
            "Sap",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh,
        dimensionedScalar(dimVelocity)
    );

    areaVectorField tau2
    (
        IOobject
        (
            "tau2",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh,
        dimensionedVector(dimPressure/dimDensity)
    );

    Info << endl;
