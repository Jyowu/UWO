/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "kEpsilon_BCS3_kande.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"

namespace Foam
{
namespace RASModels
{

template<class BasicMomentumTransportModel>
tmp<volScalarField>
kEpsilon_BCS3_kande<BasicMomentumTransportModel>::boundEpsilon()
{
    tmp<volScalarField> tCmuk2(Cmu_*sqr(k_));
    epsilon_ = max(epsilon_, tCmuk2()/(this->nutMaxCoeff_*this->nu()));
    return tCmuk2;
}


template<class BasicMomentumTransportModel>
void kEpsilon_BCS3_kande<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = boundEpsilon()/epsilon_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kEpsilon_BCS3_kande<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix(k_, dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime)
    );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kEpsilon_BCS3_kande<BasicMomentumTransportModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix(epsilon_, dimVolume*this->rho_.dimensions()*epsilon_.dimensions()/dimTime)
    );
}


template<class BasicMomentumTransportModel>
kEpsilon_BCS3_kande<BasicMomentumTransportModel>::kEpsilon_BCS3_kande
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (type, alpha, rho, U, alphaRhoPhi, phi, viscosity),
    Cmu_(dimensioned<scalar>::lookupOrAddToDict("Cmu",  this->coeffDict_, 0.09)),
    C1_ (dimensioned<scalar>::lookupOrAddToDict("C1",   this->coeffDict_, 1.44)),
    C2_ (dimensioned<scalar>::lookupOrAddToDict("C2",   this->coeffDict_, 1.92)),
    C3_ (dimensioned<scalar>::lookupOrAddToDict("C3",   this->coeffDict_, 0.0)),
    sigmak_  (dimensioned<scalar>::lookupOrAddToDict("sigmak",   this->coeffDict_, 1.0)),
    sigmaEps_(dimensioned<scalar>::lookupOrAddToDict("sigmaEps", this->coeffDict_, 1.3)),
    curvatureSwirlTools_(this->mesh_, this->coeffDict_),
    k_
    (
        IOobject(this->groupName("k"), this->runTime_.name(), this->mesh_,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE),
        this->mesh_
    ),
    epsilon_
    (
        IOobject(this->groupName("epsilon"), this->runTime_.name(), this->mesh_,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    boundEpsilon();
    if (type == typeName) this->printCoeffs(type);
}


template<class BasicMomentumTransportModel>
bool kEpsilon_BCS3_kande<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        curvatureSwirlTools_.read(this->coeffDict());
        return true;
    }
    return false;
}


template<class BasicMomentumTransportModel>
void kEpsilon_BCS3_kande<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_) return;

    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints(Foam::fvConstraints::New(this->mesh_));

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U))());

    const volTensorField gradU(fvc::grad(U));

    volScalarField::Internal G
    (
        this->GName(),
        nut()*(dev(twoSymm(gradU.v())) && gradU.v())
    );

    // Evaluate curvature/swirl correction
    volScalarField oneFrEff
    (
        IOobject("frEff", this->runTime_.name(), this->mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        this->mesh_,
        dimensionedScalar("one", dimless, 1.0)
    );

    tmp<volScalarField> tFrEff;
    const volScalarField* frEffPtr = nullptr;

    if (curvatureSwirlTools_.active())
    {
        const dimensionedScalar CmuRef("CmuRef", dimless, 0.09);
        const dimensionedScalar kSmall("kSmall", k_.dimensions(), SMALL);
        const volScalarField omegaLike
        (
            typedName("omegaLike"),
            epsilon_/(CmuRef*max(k_, kSmall))
        );

        curvatureSwirlData3 swirlData = curvatureSwirlTools_.evaluate
        (
            U, gradU, rho, k_, omegaLike, nut
        );

        tFrEff = swirlData.frEff;
        frEffPtr = &tFrEff();
    }

    const volScalarField::Internal frEffI
    (
        typedName("frEffI"),
        frEffPtr ? (*frEffPtr)() : oneFrEff()
    );

    // Corrected production: GCorr = G * frEff
    const volScalarField::Internal GCorr(typedName("GCorr"), G*frEffI);

    epsilon_.boundaryFieldRef().updateCoeffs();

    // Epsilon equation — frEff via GCorr
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha()*rho()*GCorr*epsilon_()/k_()
      - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha()*rho()*divU, epsilon_)
      - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_(), epsilon_)
      + epsilonSource()
      + fvModels.source(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvConstraints.constrain(epsilon_);
    boundEpsilon();

    // k equation — frEff via GCorr
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*GCorr
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
      + kSource()
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);

    correctNut();
}

} // End namespace RASModels
} // End namespace Foam
