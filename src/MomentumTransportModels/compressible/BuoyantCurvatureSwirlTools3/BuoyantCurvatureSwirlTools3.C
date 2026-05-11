#include "BuoyantCurvatureSwirlTools3.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMesh.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"
#include "mathematicalConstants.H"

namespace Foam
{

defineTypeNameAndDebug(BuoyantCurvatureSwirlTools3, 0);


// * * * * * * * * * * * * * * Private Helper Functions  * * * * * * * * * //

vector BuoyantCurvatureSwirlTools3::normalisedAxis(const vector& axis)
{
    scalar magAxis = mag(axis);

    if (magAxis < SMALL)
    {
        FatalErrorInFunction
            << "axisDirection has near-zero magnitude: " << axis
            << exit(FatalError);
    }

    return axis/magAxis;
}


vectorField BuoyantCurvatureSwirlTools3::cellRelPosition() const
{
    const vectorField& C = mesh_.C();

    vectorField rel(C.size(), Zero);

    forAll(C, i)
    {
        rel[i] = C[i] - axisOrigin_;
    }

    return rel;
}


tmp<volScalarField> BuoyantCurvatureSwirlTools3::axialVelocity
(
    const volVectorField& U
) const
{
    const dimensionedVector eAxis
    (
        "eAxis",
        dimless,
        axisDirection_
    );

    tmp<volScalarField> tUaxial
    (
        new volScalarField
        (
            IOobject
            (
                "Uaxial",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimVelocity, 0.0)
        )
    );

    tUaxial.ref() = (U & eAxis);

    return tUaxial;
}


tmp<volScalarField> BuoyantCurvatureSwirlTools3::tangentialVelocity
(
    const volVectorField& U
) const
{
    const vectorField relPos(cellRelPosition());

    volScalarField* uThetaPtr = new volScalarField
    (
        IOobject
        (
            "uTheta",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimVelocity, 0.0)
    );

    volScalarField& uTheta = *uThetaPtr;

    forAll(U, i)
    {
        const vector& rel = relPos[i];
        const scalar axialProj = (rel & axisDirection_);
        const vector vr = rel - axialProj*axisDirection_;
        const scalar rMag = mag(vr);
        const vector er = vr/(rMag + SMALL);
        const vector eTheta = axisDirection_ ^ er;

        uTheta[i] = (U[i] & eTheta);
    }

    uTheta.correctBoundaryConditions();

    return tmp<volScalarField>(uThetaPtr);
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

BuoyantCurvatureSwirlTools3::BuoyantCurvatureSwirlTools3
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    curvatureCorrection_
    (
        dict.lookupOrDefault<Switch>("curvatureCorrection", true)
    ),
    swirlCorrection_
    (
        dict.lookupOrDefault<Switch>("swirlCorrection", true)
    ),
    buoyancyCorrection_
    (
        dict.lookupOrDefault<Switch>("buoyancyCorrection", false)
    ),
    axisOrigin_(dict.lookupOrDefault<vector>("axisOrigin", vector::zero)),
    axisDirection_
    (
        normalisedAxis
        (
            dict.lookupOrDefault<vector>("axisDirection", vector(0, 0, 1))
        )
    ),
    frMax_
    (
        "frMax",
        dimless,
        dict.lookupOrDefault<scalar>("frMax", 1.25)
    ),
    cCurv_
    (
        "cCurv",
        dimless,
        dict.lookupOrDefault<scalar>("cCurv", 1.0)
    ),
    cr1_
    (
        "cr1",
        dimless,
        dict.lookupOrDefault<scalar>("cr1", 1.0)
    ),
    cr2_
    (
        "cr2",
        dimless,
        dict.lookupOrDefault<scalar>("cr2", 2.0)
    ),
    cr3_
    (
        "cr3",
        dimless,
        dict.lookupOrDefault<scalar>("cr3", 1.0)
    ),
    swirlEps_
    (
        "swirlEps",
        sqr(dimVelocity),
        dict.lookupOrDefault<scalar>("swirlEps", 1e-6)
    ),
    Cg_
    (
        "Cg",
        dimless,
        dict.lookupOrDefault<scalar>("Cg", 1.0)
    ),
    writeFields_(dict.lookupOrDefault<Switch>("writeFields", false))
{
    if (mag(axisDirection_) < SMALL)
    {
        FatalErrorInFunction
            << "Invalid axisDirection after normalization."
            << exit(FatalError);
    }

    if (cCurv_.value() <= 0)
    {
        FatalErrorInFunction
            << "Require cCurv > 0, but got cCurv=" << cCurv_.value()
            << exit(FatalError);
    }

    Info<< "BuoyantCurvatureSwirlTools3 (BCST3): "
        << "F_swirl = |u_theta| / sqrt(u_theta^2 + u_z^2 + eps)" << nl;
}


bool BuoyantCurvatureSwirlTools3::read(const dictionary& dict)
{
    curvatureCorrection_ =
        dict.lookupOrDefault<Switch>("curvatureCorrection", curvatureCorrection_);
    swirlCorrection_ =
        dict.lookupOrDefault<Switch>("swirlCorrection", swirlCorrection_);
    buoyancyCorrection_ =
        dict.lookupOrDefault<Switch>("buoyancyCorrection", buoyancyCorrection_);
    axisOrigin_ = dict.lookupOrDefault<vector>("axisOrigin", axisOrigin_);
    axisDirection_ = normalisedAxis
    (
        dict.lookupOrDefault<vector>("axisDirection", axisDirection_)
    );

    frMax_.readIfPresent(dict);
    cCurv_.readIfPresent(dict);
    cr1_.readIfPresent(dict);
    cr2_.readIfPresent(dict);
    cr3_.readIfPresent(dict);
    swirlEps_.readIfPresent(dict);
    Cg_.readIfPresent(dict);
    writeFields_ = dict.lookupOrDefault<Switch>("writeFields", writeFields_);

    if (mag(axisDirection_) < SMALL)
    {
        FatalErrorInFunction
            << "Invalid axisDirection after normalization."
            << exit(FatalError);
    }

    if (cCurv_.value() <= 0)
    {
        FatalErrorInFunction
            << "Require cCurv > 0, but got cCurv=" << cCurv_.value()
            << exit(FatalError);
    }

    return true;
}


bool BuoyantCurvatureSwirlTools3::active() const
{
    return curvatureCorrection_ || swirlCorrection_ || buoyancyCorrection_;
}


curvatureSwirlData3 BuoyantCurvatureSwirlTools3::evaluate
(
    const volVectorField& U,
    const volTensorField& gradU,
    const volScalarField& rho,
    const volScalarField& k,
    const volScalarField& omegaLike,
    const volScalarField& nut
) const
{
    curvatureSwirlData3 data;

    const dimensionedScalar zero(dimless, 0.0);
    const dimensionedScalar one(dimless, 1.0);
    const dimensionedScalar zeroGb("zeroGb", k.dimensions()/dimTime, 0.0);
    const bool needSwirlFields =
        swirlCorrection_ || writeFields_;
    const bool needCurvatureFields =
        curvatureCorrection_ || writeFields_;

    if (needSwirlFields)
    {
        data.uTheta = tangentialVelocity(U);
        data.Uaxial = axialVelocity(U);

        // F_swirl = |u_θ| / sqrt(u_θ² + u_z² + ε)
        data.Fswirl = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Fswirl",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                min
                (
                    max
                    (
                        mag(data.uTheta())
                      / sqrt
                        (
                            sqr(data.uTheta())
                          + sqr(data.Uaxial())
                          + swirlEps_
                        ),
                        zero
                    ),
                    one
                )
            )
        );

        if (!swirlCorrection_)
        {
            data.Fswirl.ref() = zero;
        }
    }
    else
    {
        data.Fswirl = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Fswirl",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                zero
            )
        );
    }

    if (needCurvatureFields && curvatureCorrection_)
    {
        const volSymmTensorField S(symm(gradU));
        const volTensorField Omega(skew(gradU));

        const volScalarField strainMag
        (
            typedName("strainMag"),
            sqrt(2.0*magSqr(S))
        );

        const volScalarField rotationMag
        (
            typedName("rotationMag"),
            sqrt(2.0*magSqr(Omega))
        );

        const scalar deltaTValue = max(mesh_.time().deltaTValue(), SMALL);
        const dimensionedScalar deltaT("deltaT", dimTime, deltaTValue);
        const dimensionedScalar omegaSmall
        (
            "omegaSmall",
            omegaLike.dimensions(),
            SMALL
        );
        const dimensionedScalar d2Small(sqr(omegaSmall));
        const dimensionedScalar d4Small(sqr(d2Small));

        const volScalarField D2
        (
            typedName("D2"),
            max(sqr(strainMag), 0.09*sqr(omegaLike))
        );

        const volScalarField D
        (
            typedName("D"),
            sqrt(D2)
        );

        const volSymmTensorField Sold(symm(fvc::grad(U.oldTime())));

        const volScalarField DSxxDt
        (
            typedName("DSxxDt"),
            (S.component(symmTensor::XX) - Sold.component(symmTensor::XX))/deltaT
          + (U & fvc::grad(S.component(symmTensor::XX)))
        );

        const volScalarField DSxyDt
        (
            typedName("DSxyDt"),
            (S.component(symmTensor::XY) - Sold.component(symmTensor::XY))/deltaT
          + (U & fvc::grad(S.component(symmTensor::XY)))
        );

        const volScalarField DSxzDt
        (
            typedName("DSxzDt"),
            (S.component(symmTensor::XZ) - Sold.component(symmTensor::XZ))/deltaT
          + (U & fvc::grad(S.component(symmTensor::XZ)))
        );

        const volScalarField DSyyDt
        (
            typedName("DSyyDt"),
            (S.component(symmTensor::YY) - Sold.component(symmTensor::YY))/deltaT
          + (U & fvc::grad(S.component(symmTensor::YY)))
        );

        const volScalarField DSyzDt
        (
            typedName("DSyzDt"),
            (S.component(symmTensor::YZ) - Sold.component(symmTensor::YZ))/deltaT
          + (U & fvc::grad(S.component(symmTensor::YZ)))
        );

        const volScalarField DSzzDt
        (
            typedName("DSzzDt"),
            (S.component(symmTensor::ZZ) - Sold.component(symmTensor::ZZ))/deltaT
          + (U & fvc::grad(S.component(symmTensor::ZZ)))
        );

        volSymmTensorField DSDt
        (
            IOobject
            (
                "DSDt",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedSymmTensor("zero", DSxxDt.dimensions(), symmTensor::zero)
        );

        DSDt.replace(symmTensor::XX, DSxxDt);
        DSDt.replace(symmTensor::XY, DSxyDt);
        DSDt.replace(symmTensor::XZ, DSxzDt);
        DSDt.replace(symmTensor::YY, DSyyDt);
        DSDt.replace(symmTensor::YZ, DSyzDt);
        DSDt.replace(symmTensor::ZZ, DSzzDt);

        const volScalarField rStar
        (
            typedName("rStar"),
            strainMag/max(rotationMag, omegaSmall)
        );

        const volScalarField rTilde
        (
            typedName("rTilde"),
            (twoSymm(Omega & S) && DSDt)/max(rotationMag*D2*D, d4Small)
        );

        const volScalarField fRot
        (
            typedName("fRot"),
            ((one + cr1_)*((2.0*rStar)/(one + rStar))*(one - cr3_*atan(cr2_*rTilde)))
          - cr1_
        );

        const volScalarField frTilda
        (
            typedName("frTilda"),
            max(min(fRot, frMax_), zero)
        );

        data.fr = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "fr",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                one
            )
        );

        data.fr.ref() = max(zero, one + cCurv_*(frTilda - one));
    }
    else
    {
        data.fr = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "fr",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                one
            )
        );
    }

    if (curvatureCorrection_)
    {
        if (swirlCorrection_)
        {
            data.frEff = tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "frEff",
                        mesh_.time().name(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    one + data.Fswirl()*(data.fr() - one)
                )
            );
        }
        else
        {
            data.frEff = tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "frEff",
                        mesh_.time().name(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    data.fr()
                )
            );
        }
    }
    else
    {
        data.frEff = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "frEff",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                one
            )
        );
    }

    data.Gb = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Gb",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            zeroGb
        )
    );

    data.buoyLimiter = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "buoyLimiter",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            zero
        )
    );

    if (buoyancyCorrection_)
    {
        const uniformDimensionedVectorField& g =
            mesh_.objectRegistry::lookupObject<uniformDimensionedVectorField>("g");

        if (mag(g.value()) > small)
        {
            const vector gHat(g.value()/mag(g.value()));
            const volVectorField gradRho(fvc::grad(rho));
            const dimensionedScalar rhoSmall("rhoSmall", rho.dimensions(), SMALL);
            const dimensionedScalar uSmall("uSmall", dimVelocity, SMALL);

            data.Gb.ref() =
                Cg_*nut*(g & gradRho)/max(rho, rhoSmall);

            const volScalarField v(typedName("v"), gHat & U);
            const volScalarField u
            (
                typedName("u"),
                mag(U - gHat*v) + uSmall
            );

            data.buoyLimiter.ref() = tanh(mag(v)/u);
        }
    }

    if (writeFields_ && mesh_.time().writeTime())
    {
        if (data.uTheta.valid()) data.uTheta().write();
        if (data.Uaxial.valid()) data.Uaxial().write();
        if (data.Fswirl.valid()) data.Fswirl().write();
        if (data.fr.valid()) data.fr().write();
        data.frEff().write();
        data.Gb().write();
        data.buoyLimiter().write();
    }

    return data;
}


} // End namespace Foam
