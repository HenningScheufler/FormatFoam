autoPtr<suspensionFrictionModel> friction2
(
    suspensionFrictionModel::New(transportProperties2, UsDiff, h2, c2)
);

autoPtr<suspensionEntrainmentModel> entrainment2
(
    suspensionEntrainmentModel::New(transportProperties2, Us2, h2, hentrain, c2, tau2)
);

autoPtr<ambientEntrainmentModel> ambientEntrainment2
(
    ambientEntrainmentModel::New(transportProperties2, Us2, h2, c2)
);

autoPtr<suspensionDepositionModel> sedimentation2
(
    suspensionDepositionModel::New(transportProperties2, Us2, h2, c2, tau2)
);

autoPtr<couplingModel> coupling
(
    couplingModel::New(transportProperties, Us1, h1, pb1, Us2, h2, c2)
);
