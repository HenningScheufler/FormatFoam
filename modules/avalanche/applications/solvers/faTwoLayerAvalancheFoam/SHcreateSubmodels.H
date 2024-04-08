autoPtr<frictionModel> friction1
(
    frictionModel::New(transportProperties1, Us1, h1, pb1)
);

autoPtr<entrainmentModel> entrainment1
(
    entrainmentModel::New(transportProperties1, Us1, h1, hentrain, pb1, tau1)
);

autoPtr<depositionModel> deposition1
(
    depositionModel::New(transportProperties1, Us1, h1, hentrain, pb1, tau1)
);
