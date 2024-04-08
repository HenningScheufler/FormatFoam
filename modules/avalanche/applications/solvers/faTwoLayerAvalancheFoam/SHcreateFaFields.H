    Info << "Reading field h1" << endl;
    areaScalarField h1
    (
        IOobject
        (
            "h1",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh
    );


    Info << "Reading field Us1" << endl;
    areaVectorField Us1
    (
        IOobject
        (
            "Us1",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh
    );


    edgeScalarField phis1
    (
        IOobject
        (
            "phis1",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fac::interpolate(Us1) & aMesh.Le()
    );


    edgeScalarField phi2s1
    (
        IOobject
        (
            "phi2s1",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fac::interpolate(h1*Us1) & aMesh.Le()
    );

    const areaVectorField& n = aMesh.faceAreaNormals();

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

    areaScalarField pb1
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

    areaVectorField tau1
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
