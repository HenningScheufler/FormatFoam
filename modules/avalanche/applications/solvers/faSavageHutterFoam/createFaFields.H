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

    areaScalarField pb
    (
        IOobject
        (
            "pb",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh,
        dimensionedScalar(dimPressure)
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
        dimensionedVector(dimPressure)
    );

    //entrainment height
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
        dimensionedScalar(dimLength)
    );
