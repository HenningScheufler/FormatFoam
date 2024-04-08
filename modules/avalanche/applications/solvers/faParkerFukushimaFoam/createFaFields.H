    Info << "Reading field h" << endl;
    areaScalarField h
    (
        IOobject
        (
            "h",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh
    );


    Info << "Reading field Us" << endl;
    areaVectorField Us
    (
        IOobject
        (
            "Us",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh
    );

    Info << "Reading field c" << endl;
    areaScalarField c
    (
        IOobject
        (
            "c",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh
    );

    Info << "Reading field hentrain" << endl;
    areaScalarField hentrain
    (
        IOobject
        (
            "hentrain",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        aMesh,
        dimensionedScalar("hentrain0", dimLength, GREAT)
    );

    edgeScalarField phis
    (
        IOobject
        (
            "phis",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fac::interpolate(Us) & aMesh.Le()
    );

    edgeScalarField phi2s
    (
        IOobject
        (
            "phi2s",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fac::interpolate(h*Us) & aMesh.Le()
    );


    areaVectorField n = aMesh.faceAreaNormals();

    areaScalarField gn
    (
        IOobject
        (
            "gn",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        g & n
    );

    areaVectorField gs
    (
        IOobject
        (
            "gs",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        g - gn*n
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

    areaScalarField Sd
    (
        IOobject
        (
            "Sd",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh,
        dimensionedScalar(dimVelocity)
    );

    areaScalarField Sa
    (
        IOobject
        (
            "Sa",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh,
        dimensionedScalar(dimVelocity)
    );


    areaVectorField tau
    (
        IOobject
        (
            "tau",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh,
        dimensionedVector(dimPressure/dimDensity)
    );

    areaScalarField boundaryCell
    (
        IOobject
        (
            "dist",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh,
        dimensionedScalar("one", dimless, 1)
    );

