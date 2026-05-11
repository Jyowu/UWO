/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kOmegaSST_BCS3.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"

namespace Foam
{
namespace RASModels
{

template<class BasicMomentumTransportModel>
kOmegaSST_BCS3<BasicMomentumTransportModel>::kOmegaSST_BCS3
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
    kOmegaSST<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity,
        type
    ),
    curvatureSwirlTools_(this->mesh_, this->coeffDict_)
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


template<class BasicMomentumTransportModel>
bool kOmegaSST_BCS3<BasicMomentumTransportModel>::read()
{
    if (kOmegaSST<BasicMomentumTransportModel>::read())
    {
        curvatureSwirlTools_.read(this->coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void kOmegaSST_BCS3<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    volScalarField& k = this->k_;
    volScalarField& omega = this->omega_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    kOmegaSST<BasicMomentumTransportModel>::eddyViscosity::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()()
    );

    const volTensorField gradU(fvc::grad(U));
    const volScalarField S2(typedName("S2"), 2*magSqr(symm(gradU)));
    const volScalarField::Internal GbyNu
    (
        typedName("GbyNu"),
        dev(twoSymm(gradU.v())) && gradU.v()
    );
    const volScalarField::Internal G
    (
        this->GName(),
        nut()*GbyNu
    );

    tmp<volScalarField> tFrEff;
    const volScalarField* frEffPtr = nullptr;
    tmp<volScalarField> tGbCoef;
    const volScalarField* GbCoefPtr = nullptr;
    tmp<volScalarField> tBuoyLimiter;
    const volScalarField* buoyLimiterPtr = nullptr;

    volScalarField oneFrEff
    (
        IOobject
        (
            "frEff",
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("oneFrEff", dimless, scalar(1.0))
    );

    volScalarField zeroGbCoef
    (
        IOobject
        (
            "GbCoef",
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zeroGbCoef", dimless/dimTime, 0.0)
    );

    volScalarField zeroBuoyLimiter
    (
        IOobject
        (
            "buoyLimiter",
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zeroBuoyLimiter", dimless, 0.0)
    );

    if (curvatureSwirlTools_.active())
    {
        const dimensionedScalar kSmall("kSmall", k.dimensions(), SMALL);

        curvatureSwirlData3 swirlData = curvatureSwirlTools_.evaluate
        (
            U,
            gradU,
            rho,
            k,
            omega,
            nut
        );

        tFrEff = swirlData.frEff;
        frEffPtr = &tFrEff();

        tGbCoef = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "GbCoef",
                    this->runTime_.name(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                swirlData.Gb()/(k + kSmall)
            )
        );
        GbCoefPtr = &tGbCoef();

        tBuoyLimiter = swirlData.buoyLimiter;
        buoyLimiterPtr = &tBuoyLimiter();
    }

    const volScalarField::Internal frEffI
    (
        typedName("frEffI"),
        frEffPtr ? (*frEffPtr)() : oneFrEff()
    );

    const volScalarField::Internal GCorr
    (
        this->GName(),
        G*frEffI
    );

    const volScalarField::Internal GbyNuCorr
    (
        typedName("GbyNuCorr"),
        GbyNu*frEffI
    );

    omega.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(k) & fvc::grad(omega))/omega
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega)
          + fvm::div(alphaRhoPhi, omega)
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), omega)
         ==
            alpha()*rho()*gamma
           *min
            (
                GbyNuCorr,
                (this->c1_/this->a1_)*this->betaStar_*omega()
               *max(this->a1_*omega(), this->b1_*F23()*sqrt(S2()))
            )
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega)
          - fvm::Sp(alpha()*rho()*beta*omega(), omega)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega(),
                omega
            )
          - fvm::SuSp
            (
                alpha()*rho()*gamma
               *(buoyLimiterPtr ? *buoyLimiterPtr : zeroBuoyLimiter)
               *(GbCoefPtr ? *GbCoefPtr : zeroGbCoef),
                omega
            )
          + this->Qsas(S2(), gamma, beta)
          + this->omegaSource()
          + fvModels.source(alpha, rho, omega)
        );

        omegaEqn.ref().relax();
        fvConstraints.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega.boundaryFieldRef());
        solve(omegaEqn);
        fvConstraints.constrain(omega);
        this->boundOmega();
    }

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k)
      + fvm::div(alphaRhoPhi, k)
      - fvm::laplacian(alpha*rho*this->DkEff(F1), k)
     ==
        alpha()*rho()*this->Pk(GCorr)
      - fvm::SuSp(alpha()*rho()*(GbCoefPtr ? *GbCoefPtr : zeroGbCoef), k)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k)
      - fvm::Sp(alpha()*rho()*this->epsilonByk(F1(), F23()), k)
      + this->kSource()
      + fvModels.source(alpha, rho, k)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k);
    bound(k, this->kMin_);
    this->boundOmega();

    this->correctNut(S2, F23);
}


} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
