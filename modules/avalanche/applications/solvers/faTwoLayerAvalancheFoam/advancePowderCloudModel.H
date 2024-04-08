
UsDiff = Us2-Us1;

const areaVectorField & tauSc2 = friction2->tauSc();
const areaScalarField & tauSp2 = friction2->tauSp();

geff = max(gn - (fac::ndiv(phis2, Us2)&n), dimensionedScalar(dimVelocity/dimTime));

//Parker et al. (1986), Eq. (5)
faVectorMatrix Us2Eqn
(
    fam::ddt(h2, Us2)
  + xi2*fam::div(phi2s2, Us2)
  + tauSc2                // ustar**2 expl. part
  - tauSp2*Us1            // ustar**2 expl. part
  + fam::Sp(tauSp2, Us2)  // ustar**2 impl. part
  ==
    R2*gs*c2*h2
  - 1./2.*R2*fac::grad(geff*c2*h2*h2)
  + xit1*Sdp*Us1*(1.+R2)
  - fam::Sp
    (
        Spd,
        Us2
    )
);


if (correctMomentum2)
{
    Us2Eqn += R2*c2*fam::ddt(h2, Us2)
            + R2*h2*Us2*fac::ddt(c2)
            + xi2*R2*fam::div(fac::interpolate(c2)*phi2s2, Us2);
}

if (!final)
    Us2Eqn.relax();

SolverPerformance<vector> Us2Residual = solve(Us2Eqn);

Us2.correctBoundaryConditions();

phis2 = (fac::interpolate(Us2) & aMesh.Le());
tau2 = tauSc2 + tauSp2*Us2;

//Calculate sediment entrainment
areaScalarField Sm2 = entrainment2->Sm();

//no entrainment if there is not a minimal flow thickness and not if the dense flow layer is below
Sm2 = Sm2*pos(h2-hentmin2)*powderLayerAlone;
Sm2 = min(Sm2, hentrain/runTime.deltaT());

//Calculate ambient fluid entrainment
Sap = ambientEntrainment2->Sm();
Sap = Sap*pos(hwem2-h2); //clip entrainment if maximum water height is reached

//Additional sink term for water
const areaScalarField vs(R2*geff*Ds2*Ds2/18./nu2); //Terminal velocity
const areaScalarField dw(pos(h2-hmin2*2.)*vs);


//Parker et al. (1986), Eq. (3)
faScalarMatrix h2Eqn
(
    fam::ddt(h2)
  + fam::div(phis2, h2)
  ==
    Sap
  + Sdp/(h1 + dimensionedScalar("small", dimLength, 1e-5)) * h1 / packingd2
  - fam::Sp
    (
        Spd/(h2 + dimensionedScalar("small", dimLength, 1e-5)) / packingd2,
        h2
    )
);

if (waterSink2)
{
    h2Eqn -= fam::Sp(dw/(h2+dimensionedScalar("hmin", dimLength, 1e-5)),h2);
}

solverPerformance h2Residual = h2Eqn.solve();

phi2s2 = h2Eqn.flux();

if (!final)
{
     h2.relax();
}

if (bindHeight2)
{
    h2 = max(h2, hmin2);
}


faScalarMatrix c2Eqn
(
    fam::ddt(h2, c2)
  + fam::div(phi2s2, c2)
 ==
    Sm2/(hentrain+dimensionedScalar("hentrainmin", dimLength, 1e-5))*hentrain
  + Sdp/(h1 + dimensionedScalar("small", dimLength, 1e-5)) * h1
  - fam::Sp(Spd/(c2+dimensionedScalar("cmin", dimless, 1e-7)), c2)
);

if(!final)
    c2Eqn.relax();

solverPerformance c2Residual = c2Eqn.solve();

if (bindHeight2)
    c2 = max(c2, cmin2);

