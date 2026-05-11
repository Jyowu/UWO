#ifndef BuoyantCurvatureSwirlTools3_H
#define BuoyantCurvatureSwirlTools3_H

#include "volFields.H"
#include "fvMesh.H"
#include "dictionary.H"
#include "Switch.H"
#include "vector.H"
#include "tmp.H"

namespace Foam
{

struct curvatureSwirlData3
{
    tmp<volScalarField> uTheta;
    tmp<volScalarField> Uaxial;
    tmp<volScalarField> Fswirl;
    tmp<volScalarField> fr;
    tmp<volScalarField> frEff;
    tmp<volScalarField> Gb;
    tmp<volScalarField> buoyLimiter;
};


class BuoyantCurvatureSwirlTools3
{
    // Private data

        const fvMesh& mesh_;

        Switch curvatureCorrection_;
        Switch swirlCorrection_;
        Switch buoyancyCorrection_;
        vector axisOrigin_;
        vector axisDirection_;

        dimensionedScalar frMax_;
        dimensionedScalar cCurv_;
        dimensionedScalar cr1_;
        dimensionedScalar cr2_;
        dimensionedScalar cr3_;
        dimensionedScalar swirlEps_;
        dimensionedScalar Cg_;

        Switch writeFields_;


    // Private Member Functions

        static vector normalisedAxis(const vector& axis);

        vectorField cellRelPosition() const;

        tmp<volScalarField> axialVelocity
        (
            const volVectorField& U
        ) const;

        tmp<volScalarField> tangentialVelocity
        (
            const volVectorField& U
        ) const;


public:

    TypeName("BuoyantCurvatureSwirlTools3");


    // Constructors

        BuoyantCurvatureSwirlTools3
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        virtual ~BuoyantCurvatureSwirlTools3() = default;


    // Member Functions

        bool active() const;

        curvatureSwirlData3 evaluate
        (
            const volVectorField& U,
            const volTensorField& gradU,
            const volScalarField& rho,
            const volScalarField& k,
            const volScalarField& omegaLike,
            const volScalarField& nut
        ) const;

        bool read(const dictionary& dict);
};


} // End namespace Foam

#endif
