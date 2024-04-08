
autoPtr<suspensionEntrainmentModel> entrainment
(
    suspensionEntrainmentModel::New(transportProperties, Us, h, hentrain, c, tau)
);


autoPtr<ambientEntrainmentModel> ambientEntrainment
(
    ambientEntrainmentModel::New(transportProperties, Us, h, c)
);


autoPtr<suspensionFrictionModel> friction
(
    suspensionFrictionModel::New(transportProperties, Us, h, c)
);


autoPtr<suspensionDepositionModel> sedimentation
(
    suspensionDepositionModel::New(transportProperties, Us, h, c, tau)
);


if (relaxBoundaries)
{
    const faPatchList& patches = aMesh.boundary();

    const labelList& owner = aMesh.edgeOwner();

    forAll(Us.boundaryField(), patchi)
    {
        if (Us.boundaryField()[patchi].type() != "processor" &&
            Us.boundaryField()[patchi].fixesValue()){

            for(label faceI = patches[patchi].start(); faceI < patches[patchi].start()+patches[patchi].size(); faceI++)
            {
                boundaryCell[owner[faceI]] = 0;
            }
        }
    }
}
