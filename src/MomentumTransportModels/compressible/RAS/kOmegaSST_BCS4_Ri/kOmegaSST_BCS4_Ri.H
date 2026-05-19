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

Class
    Foam::RASModels::kOmegaSST_BCS4_Ri

Description
    k-omega SST turbulence model using BuoyantCurvatureSwirlTools4_Ri for
    curvature, swirl and optional buoyancy modifications.

    Curvature correction uses the simplified Richardson-number formulation
    (Alahmadi et al. 2021):
        r~ = (Omega/S)(Omega/S - 1)
        r* = S/Omega
    avoiding Lagrangian derivatives of the strain-rate tensor.

    Swirl blending:
        F_swirl = |u_theta| / sqrt(u_theta^2 + u_z^2 + eps)

SourceFiles
    kOmegaSST_BCS4_Ri.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSST_BCS4_Ri_H
#define kOmegaSST_BCS4_Ri_H

#include "kOmegaSST.H"
#include "BuoyantCurvatureSwirlTools4_Ri.H"

namespace Foam
{
namespace RASModels
{

template<class BasicMomentumTransportModel>
class kOmegaSST_BCS4_Ri
:
    public kOmegaSST<BasicMomentumTransportModel>
{
protected:

    BuoyantCurvatureSwirlTools4_Ri curvatureSwirlTools_;

    volScalarField G_;

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;

    TypeName("kOmegaSST_BCS4_Ri");

    kOmegaSST_BCS4_Ri
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const viscosity& viscosity,
        const word& type = typeName
    );

    virtual ~kOmegaSST_BCS4_Ri()
    {}

    virtual bool read();

    virtual void correct();
};


} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "kOmegaSST_BCS4_Ri.C"
#endif

#endif

// ************************************************************************* //
