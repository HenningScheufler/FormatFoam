
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

dictionary transportProperties1 = transportProperties.subDict("denseFlow");

dimensionedScalar hmin1("hmin", dimLength, transportProperties1);

dimensionedScalar xi1("xi", dimless, transportProperties1);

dimensionedScalar xit1("xit", dimless, transportProperties1);

dimensionedScalar packingd("phi", dimless, transportProperties1);

Switch pressureFeedback1(transportProperties1.get<Switch>("pressureFeedback"));

Switch bindHeight1(transportProperties1.getOrDefault<Switch>("bindHeight", true));

Info << "Running dense flow with" << endl
     << "    xi " << xi1 << endl
     << "    xi_t" << xit1 << endl
     << "    hmin " << hmin1 << endl
     << "    pressureFeedback is " << pressureFeedback1 << endl
     << "    bindHeight is " << bindHeight1 << endl << endl;
