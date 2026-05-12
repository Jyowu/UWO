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
    Foam::RASModels::kOmegaSST_BCS3

Description
    k-omega SST turbulence model using the shared
    BuoyantCurvatureSwirlTools3 helper for curvature, swirl and optional
    buoyancy modifications.

    Swirl blending uses the velocity-based formulation:
        F_swirl = |u_theta| / sqrt(u_theta^2 + u_z^2 + eps)

SourceFiles
    kOmegaSST_BCS3.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSST_BCS3_H
#define kOmegaSST_BCS3_H

#include "kOmegaSST.H"
#include "BuoyantCurvatureSwirlTools3.H"

namespace Foam
{
namespace RASModels
{

template<class BasicMomentumTransportModel>
class kOmegaSST_BCS3
:
    public kOmegaSST<BasicMomentumTransportModel>
{
protected:

    BuoyantCurvatureSwirlTools3 curvatureSwirlTools_;

    volScalarField G_;

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;

    TypeName("kOmegaSST_BCS3");

    kOmegaSST_BCS3
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const viscosity& viscosity,
        const word& type = typeName
    );

    virtual ~kOmegaSST_BCS3()
    {}

    virtual bool read();

    virtual void correct();
};


} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "kOmegaSST_BCS3.C"
#endif

#endif

// ************************************************************************* //
