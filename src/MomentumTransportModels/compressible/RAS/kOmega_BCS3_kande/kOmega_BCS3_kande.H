/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::RASModels::kOmega_BCS3_kande

Description
    Standard k-omega with BuoyantCurvatureSwirlTools3 correction applied
    to both k and omega production (k-and-e).

    frEff scales G in both equations:
        k     : alpha*rho*G*frEff
        omega : gamma*alpha*rho*G*frEff*omega/k

    Default coefficients identical to standard kOmega.

SourceFiles
    kOmega_BCS3_kande.C

\*---------------------------------------------------------------------------*/

#ifndef kOmega_BCS3_kande_H
#define kOmega_BCS3_kande_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "BuoyantCurvatureSwirlTools3.H"

namespace Foam
{
namespace RASModels
{

template<class BasicMomentumTransportModel>
class kOmega_BCS3_kande
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel>>
{
protected:

    dimensionedScalar betaStar_;
    dimensionedScalar beta_;
    dimensionedScalar gamma_;
    dimensionedScalar alphaK_;
    dimensionedScalar alphaOmega_;

    BuoyantCurvatureSwirlTools3 curvatureSwirlTools_;

    volScalarField k_;
    volScalarField omega_;

    void boundOmega();
    virtual void correctNut();
    virtual tmp<fvScalarMatrix> kSource() const;
    virtual tmp<fvScalarMatrix> omegaSource() const;

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;

    TypeName("kOmega_BCS3_kande");

    kOmega_BCS3_kande
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const viscosity& viscosity,
        const word& type = typeName
    );

    virtual ~kOmega_BCS3_kande() {}

    virtual bool read();

    tmp<volScalarField> DkEff() const
    {
        return volScalarField::New("DkEff", alphaK_*this->nut_ + this->nu());
    }

    tmp<volScalarField> DomegaEff() const
    {
        return volScalarField::New("DomegaEff", alphaOmega_*this->nut_ + this->nu());
    }

    virtual tmp<volScalarField> k() const { return k_; }

    virtual tmp<volScalarField> omega() const { return omega_; }

    virtual tmp<volScalarField> epsilon() const
    {
        return volScalarField::New("epsilon", betaStar_*k_*omega_);
    }

    virtual void correct();
};

} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "kOmega_BCS3_kande.C"
#endif

#endif
