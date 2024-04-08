dictionary transportProperties2 = transportProperties.subDict("powderCloud");

dimensionedScalar hmin2("hmin", dimLength, transportProperties2);

dimensionedScalar cmin2("cmin", dimless, transportProperties2);

dimensionedScalar h02("h0", dimLength, transportProperties2);

dimensionedScalar u02("u0", dimVelocity, transportProperties2);

dimensionedScalar xi2("xi", dimless, transportProperties2);

dimensionedScalar R2("R", dimless, transportProperties2);

dimensionedScalar Ds2("Ds", dimLength, transportProperties2);

dimensionedScalar nu2("nu", dimViscosity, transportProperties2);

dimensionedScalar packingd2("phi0", dimless, transportProperties2);

dimensionedScalar hwem2("hwem", dimLength, transportProperties2);

dimensionedScalar hentmin2("hentmin", dimLength, transportProperties2);

Switch waterSink2(transportProperties2.getOrDefault<Switch>("waterSink", false));

Switch bindHeight2(transportProperties.getOrDefault<Switch>("bindHeight", true));

Switch correctMomentum2(transportProperties.getOrDefault<Switch>("correctMomentum", true));

Info << "Running powder cloud with" << endl
     << "    xi2 " << xi1 << endl
     << "    hwem " << hwem2 << endl
     << "    hmin " << hmin2 << endl
     << "    cmin " << cmin2 << endl
     << "    h0   " << h02 << endl
     << "    u0   " << u02 << endl
     << "    waterSink is " << waterSink2 << endl
     << "    bindHeight is " << bindHeight2 << endl << endl;

