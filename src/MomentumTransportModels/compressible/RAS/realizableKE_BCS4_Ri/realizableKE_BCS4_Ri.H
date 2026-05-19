/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

Class
    Foam::RASModels::realizableKE_BCS4_Ri

Description
    Realizable k-epsilon turbulence model using BuoyantCurvatureSwirlTools4_Ri
    for curvature, swirl and optional buoyancy modifications.

    Curvature correction uses the simplified Richardson-number formulation
    (Alahmadi et al. 2021):
        r~ = (Omega/S)(Omega/S - 1)
        r* = S/Omega
    avoiding Lagrangian derivatives of the strain-rate tensor.

    Swirl blending:
        F_swirl = |u_theta| / sqrt(u_theta^2 + u_z^2 + eps)

    frEff is applied to the k production term only (epsilon equation uses
    C1*S*epsilon form which does not carry a direct G term).

    The default model coefficients are
    \verbatim
        realizableKE_BCS4_RiCoeffs
        {
            A0          4.0;
            C2          1.9;
            sigmak      1.0;
            sigmaEps    1.2;
        }
    \endverbatim

SourceFiles
    realizableKE_BCS4_Ri.C

\*---------------------------------------------------------------------------*/

#ifndef realizableKE_BCS4_Ri_H
#define realizableKE_BCS4_Ri_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "BuoyantCurvatureSwirlTools4_Ri.H"

namespace Foam
{
namespace RASModels
{

template<class BasicMomentumTransportModel>
class realizableKE_BCS4_Ri
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel>>
{
protected:

    dimensionedScalar A0_;
    dimensionedScalar C2_;
    dimensionedScalar sigmak_;
    dimensionedScalar sigmaEps_;

    BuoyantCurvatureSwirlTools4_Ri curvatureSwirlTools_;

    volScalarField k_;
    volScalarField epsilon_;
    volScalarField G_;

    void boundEpsilon();

    tmp<volScalarField> rCmu
    (
        const volTensorField& gradU,
        const volScalarField& S2,
        const volScalarField& magS
    );

    virtual void correctNut
    (
        const volTensorField& gradU,
        const volScalarField& S2,
        const volScalarField& magS
    );

    virtual void correctNut();

    virtual tmp<fvScalarMatrix> kSource() const;

    virtual tmp<fvScalarMatrix> epsilonSource() const;

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;

    TypeName("realizableKE_BCS4_Ri");

    realizableKE_BCS4_Ri
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const viscosity& viscosity,
        const word& type = typeName
    );

    virtual ~realizableKE_BCS4_Ri()
    {}

    virtual bool read();

    tmp<volScalarField> DkEff() const
    {
        return volScalarField::New
        (
            "DkEff",
            (this->nut_/sigmak_ + this->nu())
        );
    }

    tmp<volScalarField> DepsilonEff() const
    {
        return volScalarField::New
        (
            "DepsilonEff",
            (this->nut_/sigmaEps_ + this->nu())
        );
    }

    virtual tmp<volScalarField> k() const
    {
        return k_;
    }

    virtual tmp<volScalarField> epsilon() const
    {
        return epsilon_;
    }

    virtual tmp<volScalarField> omega() const
    {
        return volScalarField::New
        (
            "omega",
            epsilon_/(0.09*k_)
        );
    }

    virtual void correct();
};


} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "realizableKE_BCS4_Ri.C"
#endif

#endif

// ************************************************************************* //
