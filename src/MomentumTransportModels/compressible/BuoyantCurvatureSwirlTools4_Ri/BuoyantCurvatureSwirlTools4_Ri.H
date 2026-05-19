#ifndef BuoyantCurvatureSwirlTools4_Ri_H
#define BuoyantCurvatureSwirlTools4_Ri_H

#include "volFields.H"
#include "fvMesh.H"
#include "dictionary.H"
#include "Switch.H"
#include "vector.H"
#include "tmp.H"

namespace Foam
{

struct curvatureSwirlData4Ri
{
    tmp<volScalarField> uTheta;
    tmp<volScalarField> Uaxial;
    tmp<volScalarField> Fswirl;
    tmp<volScalarField> fr;
    tmp<volScalarField> frEff;
    tmp<volScalarField> Gb;
    tmp<volScalarField> buoyLimiter;
};


class BuoyantCurvatureSwirlTools4_Ri
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

    TypeName("BuoyantCurvatureSwirlTools4_Ri");


    // Constructors

        BuoyantCurvatureSwirlTools4_Ri
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        virtual ~BuoyantCurvatureSwirlTools4_Ri() = default;


    // Member Functions

        bool active() const;

        curvatureSwirlData4Ri evaluate
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
