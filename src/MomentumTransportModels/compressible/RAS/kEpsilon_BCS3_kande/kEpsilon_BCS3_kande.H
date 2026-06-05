/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::RASModels::kEpsilon_BCS3_kande

Description
    Standard k-epsilon with BuoyantCurvatureSwirlTools3 correction applied
    to both k and epsilon production (k-and-e).

    frEff scales G in both equations:
        k   : alpha*rho*G*frEff
        eps : C1*alpha*rho*G*frEff*epsilon/k

    Default coefficients identical to standard kEpsilon.

SourceFiles
    kEpsilon_BCS3_kande.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilon_BCS3_kande_H
#define kEpsilon_BCS3_kande_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "BuoyantCurvatureSwirlTools3.H"

namespace Foam
{
namespace RASModels
{

template<class BasicMomentumTransportModel>
class kEpsilon_BCS3_kande
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel>>
{
protected:

    dimensionedScalar Cmu_;
    dimensionedScalar C1_;
    dimensionedScalar C2_;
    dimensionedScalar C3_;
    dimensionedScalar sigmak_;
    dimensionedScalar sigmaEps_;

    BuoyantCurvatureSwirlTools3 curvatureSwirlTools_;

    volScalarField k_;
    volScalarField epsilon_;

    tmp<volScalarField> boundEpsilon();
    virtual void correctNut();
    virtual tmp<fvScalarMatrix> kSource() const;
    virtual tmp<fvScalarMatrix> epsilonSource() const;

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;

    TypeName("kEpsilon_BCS3_kande");

    kEpsilon_BCS3_kande
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const viscosity& viscosity,
        const word& type = typeName
    );

    virtual ~kEpsilon_BCS3_kande() {}

    virtual bool read();

    tmp<volScalarField> DkEff() const
    {
        return volScalarField::New("DkEff", this->nut_/sigmak_ + this->nu());
    }

    tmp<volScalarField> DepsilonEff() const
    {
        return volScalarField::New("DepsilonEff", this->nut_/sigmaEps_ + this->nu());
    }

    virtual tmp<volScalarField> k() const { return k_; }

    virtual tmp<volScalarField> epsilon() const { return epsilon_; }

    virtual tmp<volScalarField> omega() const
    {
        return volScalarField::New("omega", epsilon_/(Cmu_*k_));
    }

    virtual void correct();
};

} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "kEpsilon_BCS3_kande.C"
#endif

#endif
