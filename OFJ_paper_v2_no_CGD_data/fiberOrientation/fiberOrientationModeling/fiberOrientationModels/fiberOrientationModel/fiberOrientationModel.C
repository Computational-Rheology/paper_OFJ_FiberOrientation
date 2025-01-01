/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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

#include "fiberOrientationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fiberOrientation
    {
        defineTypeNameAndDebug(fiberOrientationModel, 0);

        defineRunTimeSelectionTable
        (
            fiberOrientationModel, // Current class name
            dictionary             // dictionary
        );
    }
}

// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fiberOrientation::fiberOrientationModel::computeXi() const
{
	const scalar ar = readScalar(fiberPropertiesDict_.lookup("ar"));
	return (ar*ar - 1.0)/(ar*ar + 1.0);
}


Foam::scalar Foam::fiberOrientation::fiberOrientationModel::fiberVolumeFraction() const
{
    dimensionedScalar VFraction("VFraction", dimless, Zero);

    if (fiberPropertiesDict_.found("volumeFraction"))
    {
        VFraction.value() = readScalar(fiberPropertiesDict_.lookup("volumeFraction"));
    }
    else
    {
        const dimensionedScalar polymerDensity  ("polymerDensity", dimDensity, fiberPropertiesDict_);
        const dimensionedScalar fiberDensity    ("fiberDensity", dimDensity, fiberPropertiesDict_);
        const dimensionedScalar fiberWeightFrac ("fiberWeightFraction", dimless, fiberPropertiesDict_);

        VFraction = (fiberWeightFrac/fiberDensity)/(
                                                        (fiberWeightFrac/fiberDensity) + ((1.0 - fiberWeightFrac)/polymerDensity)
                                                    );
    }

    return VFraction.value();
}


Foam::volSymmTensorField& Foam::fiberOrientation::fiberOrientationModel::createA2()
{
    const word fieldName ("A2");

    if (!mesh_.foundObject<volSymmTensorField>(fieldName))
    {
        auto tfldPtr = tmp<volSymmTensorField>::New
        (
            IOobject
            (
                fieldName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        );

        mesh_.objectRegistry::store(tfldPtr);
    }

    return mesh_.lookupObjectRef<volSymmTensorField>(fieldName);
}


const Foam::volScalarField& 
Foam::fiberOrientation::fiberOrientationModel::constLookupOrConstruct(const word& name)
{
    if (!mesh_.foundObject<volScalarField>(name))
    {
        auto tfldPtr = tmp<volScalarField>::New
        (
            IOobject
            (
                name,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, One)
        );

        mesh_.objectRegistry::store(tfldPtr);
    }

    return mesh_.lookupObject<volScalarField>(name);
}


// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fiberOrientation::fiberOrientationModel::computeCI() const
{
    const scalar volumeFraction = fiberVolumeFraction();
    const scalar ar = readScalar(fiberPropertiesDict_.lookup("ar"));
    scalar CI = 0.0;

	if (volumeFraction*ar > 1.3)
    {
        // B.S. Bay, Fibre Orientation in Injection Molded Composites: A Comparison of Theory and Experiment, PhD Thesis, Mechanical Engineering, University of Illinois at Urbana-Champaign, 1991.
		CI = 0.0184*Foam::exp(-0.7148*ar*volumeFraction);
	}
	else
    {
	    // N. Phan-Thien, X. J. Fan, R. I. Tanner, and R. Zheng, “Folgar-Tucker constant for a fibre suspension in a Newtonian fluid,” J. Nonnewton. Fluid Mech., vol. 103, no. 2–3, pp. 251–260, 2002, doi: 10.1016/S0377-0257(02)00006-X.
		CI = 0.03*(1.0-Foam::exp(-0.224*ar*volumeFraction));
	}

    return CI;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fiberOrientation::fiberOrientationModel::fiberOrientationModel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    mesh_(mesh),
    fiberPropertiesDict_(dict.subDict("fiberOrientationProperties")),
    U_(U),
    phi_(phi),
    alpha_
    ( 
        constLookupOrConstruct(dict.getOrDefault<word>("phase", "none"))
    ),
    A2_(createA2()), 
    L_
    (
        IOobject
        (
            "fiberModel_L",
			mesh.time().timeName(),
			mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false    // This is needed to not have a memory leak
        ),
        fvc::grad(U)
    ),
    D_
    (
        IOobject
        (
            "fiberModel_D",
			mesh.time().timeName(),
			mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        symm(L_) 
    ),
    W_
    (
        IOobject
        (
            "fiberModel_W",
			mesh.time().timeName(),
			mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        -1.0*skew(L_) 
    ),
    shrRate_
    (
        IOobject
        (
            "fiberModel_shrRate",
			mesh.time().timeName(),
			mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        Foam::sqrt(2.0*(D_ && D_))
    ),
    D_doubleDot_A4_ 
	(
		IOobject
		(
			"fiberModel_D_doubleDot_A4",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
		),
		mesh,
        dimensionedSymmTensor("D_doubleDot_A4_", dimless/dimTime, Foam::Zero) // Made to cancel the the ddt(A_ij) 
	),
    xi_(
         fiberPropertiesDict_.found("xi") ? 
         readScalar(fiberPropertiesDict_.lookup("xi")) 
         : computeXi() 
        ),
    alphaCutOff_(dict.getOrDefault<scalar>("alphaCutOff", 0.0)),
    closureModel_(fiberOrientation::closureModel::New(fiberPropertiesDict_, A2_)),
    solutionPropertiesDict_(dict.subDict("solutionProperties")),
    nCorr_ (solutionPropertiesDict_.getOrDefault<label>("nCorr", 1)),
    schemesField_(solutionPropertiesDict_.getOrDefault("schemesField", word("A2"))),
    absTol_(solutionPropertiesDict_.getOrDefault<scalar>("absTol", 1e-5)),
    normalize_(solutionPropertiesDict_.getOrDefault<Switch>("normalize", false) ),
    updateFlowFields_(solutionPropertiesDict_.getOrDefault<Switch>("updateFlowFields", true) )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fiberOrientation::fiberOrientationModel::normalizeA2()
{     
    A2_ = A2_/tr(A2_);
    A2_.correctBoundaryConditions();
}

void Foam::fiberOrientation::fiberOrientationModel::correctFlowFieldTensors()
{
    if (updateFlowFields_)
    {
        L_ = fvc::grad(U_);

        D_ = symm(L_);

        W_ = -1.0*skew(L_);

        shrRate_ = Foam::sqrt( 2.0*(D_ && D_) );
    }
}

bool Foam::fiberOrientation::fiberOrientationModel::read()
{
    solutionPropertiesDict_.readEntry("nCorr", nCorr_);
    solutionPropertiesDict_.readEntry("absTol", absTol_);
    solutionPropertiesDict_.readIfPresent("normalize", normalize_);
    solutionPropertiesDict_.readIfPresent("updateFlowFields", updateFlowFields_);

    return true;
}

const Foam::dictionary&  Foam::fiberOrientation::fiberOrientationModel::fiberDict() const
{
    return fiberPropertiesDict_;
}


void Foam::fiberOrientation::fiberOrientationModel::write(Ostream& os) const
{
    dictionary dict;
    writeWithSubDict(os, dict);
}

void Foam::fiberOrientation::fiberOrientationModel::writeWithSubDict
(
    Ostream& os, 
    const dictionary& subDict
) const
{
    dictionary dict1("solutionProperties");
    dict1.add("nCorr", nCorr_);
    dict1.add("absTol", absTol_);
    dict1.add("normalize", normalize_);
    dict1.add("updateFlowFields", updateFlowFields_);

    dictionary dict2("fiberOrientationProperties");
    if (fiberPropertiesDict_.found("ar"))
        dict2.add("ar", name(readScalar(fiberPropertiesDict_.lookup("ar"))));

    dict2.add("model", fiberPropertiesDict_.get<word>("model"));

    dict2.add("xi", xi_);
    dict2.add("closureModel", fiberPropertiesDict_.get<word>("closureModel"));

    if (fiberPropertiesDict_.found("volumeFraction"))
        dict2.add("volumeFraction", name(readScalar(fiberPropertiesDict_.lookup("volumeFraction"))));

    dict2.add(subDict.dictName(), subDict);

    os  << "FiberOrientation Parameters: " << nl << endl;
    
    os  << indent << dict1.dictName() << dict1 << nl << endl;
    os  << indent << dict2.dictName() << dict2 << nl << endl;
}


Foam::Ostream& Foam::fiberOrientation::operator<<
(
    Ostream& os,
    const fiberOrientationModel& fm
)
{
    fm.write(os);
    os.check("Ostream& operator<<(Ostream& os, const fiberOrientationModel& fm)");
    return os;
}

// ************************************************************************* //
 