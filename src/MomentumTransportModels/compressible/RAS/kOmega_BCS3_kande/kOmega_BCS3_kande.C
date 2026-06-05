/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "kOmega_BCS3_kande.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"

namespace Foam
{
namespace RASModels
{

template<class BasicMomentumTransportModel>
void kOmega_BCS3_kande<BasicMomentumTransportModel>::boundOmega()
{
    omega_ = max(omega_, k_/(this->nutMaxCoeff_*this->nu()));
}


template<class BasicMomentumTransportModel>
void kOmega_BCS3_kande<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = k_/omega_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kOmega_BCS3_kande<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix(k_, dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime)
    );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kOmega_BCS3_kande<BasicMomentumTransportModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix(omega_, dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime)
    );
}


template<class BasicMomentumTransportModel>
kOmega_BCS3_kande<BasicMomentumTransportModel>::kOmega_BCS3_kande
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
    betaStar_  (dimensioned<scalar>::lookupOrAddToDict("betaStar",   this->coeffDict_, 0.09)),
    beta_      (dimensioned<scalar>::lookupOrAddToDict("beta",       this->coeffDict_, 0.072)),
    gamma_     (dimensioned<scalar>::lookupOrAddToDict("gamma",      this->coeffDict_, 0.52)),
    alphaK_    (dimensioned<scalar>::lookupOrAddToDict("alphaK",     this->coeffDict_, 0.5)),
    alphaOmega_(dimensioned<scalar>::lookupOrAddToDict("alphaOmega", this->coeffDict_, 0.5)),
    curvatureSwirlTools_(this->mesh_, this->coeffDict_),
    k_
    (
        IOobject(this->groupName("k"), this->runTime_.name(), this->mesh_,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE),
        this->mesh_
    ),
    omega_
    (
        IOobject(this->groupName("omega"), this->runTime_.name(), this->mesh_,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    boundOmega();
    if (type == typeName) this->printCoeffs(type);
}


template<class BasicMomentumTransportModel>
bool kOmega_BCS3_kande<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        betaStar_.readIfPresent(this->coeffDict());
        beta_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
        curvatureSwirlTools_.read(this->coeffDict());
        return true;
    }
    return false;
}


template<class BasicMomentumTransportModel>
void kOmega_BCS3_kande<BasicMomentumTransportModel>::correct()
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

    volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U))().v());

    const volTensorField gradU(fvc::grad(U));

    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(gradU.v())) && gradU.v())
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
        // Use omega_ directly as omegaLike for k-omega
        curvatureSwirlData3 swirlData = curvatureSwirlTools_.evaluate
        (
            U, gradU, rho, k_, omega_, nut
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

    omega_.boundaryFieldRef().updateCoeffs();

    // Omega equation — frEff via GCorr
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        gamma_*alpha()*rho()*GCorr*omega_()/k_()
      - fvm::SuSp(((2.0/3.0)*gamma_)*alpha()*rho()*divU, omega_)
      - fvm::Sp(beta_*alpha()*rho()*omega_(), omega_)
      + omegaSource()
      + fvModels.source(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvConstraints.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvConstraints.constrain(omega_);
    boundOmega();

    // k equation — frEff via GCorr
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*GCorr
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(betaStar_*alpha()*rho()*omega_(), k_)
      + kSource()
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);
    boundOmega();

    correctNut();
}

} // End namespace RASModels
} // End namespace Foam
