#include "calcBasalstress.H"

const areaVectorField & tauSc1 = friction1->tauSc();
const areaScalarField & tauSp1 = friction1->tauSp();

faVectorMatrix Us1Eqn
(
    fam::ddt(h1, Us1)
  + xi1*fam::div(phi2s1, Us1)
  + tauSc1
  + fam::Sp(tauSp1, Us1)
 ==
    gs*h1
  - fac::grad(pb1*h1/(2.*friction1->rho()))
  - fam::Sp
    (
        xit1*Sdp/packingd,
        Us1
    )
  + Spd*Us2
);


if (!final)
    Us1Eqn.relax();

SolverPerformance<vector> Us1Residual = solve(Us1Eqn);

tau1 = (tauSc1 + tauSp1*Us1)*friction1->rho();
phis1 = (fac::interpolate(Us1) & aMesh.Le());

const areaScalarField & Sm1 = entrainment1->Sm();
const areaScalarField & Sd1 = deposition1->Sd();


h1.storePrevIter();

faScalarMatrix h1Eqn
(
    fam::ddt(h1)
  + fam::div(phis1, h1)
 ==
    Sm1/(hentrain+dimensionedScalar("hentrainmin", dimLength, 1e-5))*hentrain
  - fam::Sp
    (
        Sd1/(h1 + dimensionedScalar("small", dimLength, 1e-5)),
        h1
    )
  + Spd/(c2+dimensionedScalar("cmin", dimless, 1e-7))*c2/packingd
  - fam::Sp
    (
        Sdp/(h1 + dimensionedScalar("small", dimLength, 1e-5))/packingd,
        h1
    )
);

solverPerformance h1Residual = h1Eqn.solve();

phi2s1 = h1Eqn.flux();

if (!final)
    h1.relax();


// Bind h
if (bindHeight1)
{
    h1 = max
    (
        h1,
        hmin1
    );
}
h1.correctBoundaryConditions();
Us1.correctBoundaryConditions();
