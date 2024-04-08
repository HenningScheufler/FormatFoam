
pb1 = -friction1->rho()*xi1*fac::ndiv(phi2s1, Us1)&n;

if (pressureFeedback1)
{
    pb1 -= fac::ngrad(pb1*h1/2.)&n;
}

pb1 += friction1->rho()*(g*h1)&n;

// Clipping
pb1 = Foam::max(pb1, dimensionedScalar(dimPressure));

