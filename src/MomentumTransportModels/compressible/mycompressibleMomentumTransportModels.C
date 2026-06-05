/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "makeCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
makeBaseMomentumTransportModel
(
    geometricOneField,
    volScalarField,
    compressibleMomentumTransportModel
);
*/

// -------------------------------------------------------------------------- //
// Laminar models
// -------------------------------------------------------------------------- //



// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //



#include "buoyantKEpsilonCC.H"
makeRASModel(buoyantKEpsilonCC);

#include "realizableKECC.H"
makeRASModel(realizableKECC);

#include "realizableKE_BCS.H"
makeRASModel(realizableKE_BCS);

#include "realizableKE_BCS1.H"
makeRASModel(realizableKE_BCS1);

#include "realizableKE_BCS2.H"
makeRASModel(realizableKE_BCS2);

#include "kOmegaSST_BCS2.H"
makeRASModel(kOmegaSST_BCS2);

#include "realizableKE_BCS3.H"
makeRASModel(realizableKE_BCS3);

#include "kOmegaSST_BCS3.H"
makeRASModel(kOmegaSST_BCS3);

#include "realizableKE_BCS4_Ri.H"
makeRASModel(realizableKE_BCS4_Ri);

#include "kOmegaSST_BCS4_Ri.H"
makeRASModel(kOmegaSST_BCS4_Ri);

#include "realizableKE_BCS3_kande.H"
makeRASModel(realizableKE_BCS3_kande);
// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //



// ************************************************************************* //
