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

#include "MRD.H"
#include "addToRunTimeSelectionTable.H"
        
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{
    defineTypeNameAndDebug(MRD, 0);

    addToRunTimeSelectionTable                                             
    (                                                                      
        fiberOrientationModel,      // Mother Class                                        
        MRD,                        // Current Class                                     
        dictionary                 // dictionary                                           
    );  
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fiberOrientation::MRD::MRD
( 
    const dictionary& dict,
    const fvMesh& mesh,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
fiberOrientationModel(dict, mesh, U, phi),
modelProperties_(fiberDict().subDict(typeName + "Properties")),
CI_(
        modelProperties_.found("CI") 
        ? 
        readScalar(modelProperties_.lookup("CI")) 
        : 
        computeCI() 
    ),
D1_(modelProperties_.getOrDefault<scalar>("D1", 1.0)),
D2_(modelProperties_.getOrDefault<scalar>("D2", 0.8)),
D3_(modelProperties_.getOrDefault<scalar>("D3", 0.15)),
C_
(
    IOobject
    (
        "fiberModel_C",
        U.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimensionedSymmTensor(dimless, Foam::Zero)
)
{
    Info << *this << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fiberOrientation::MRD::solve()
{   
    correctFlowFieldTensors();

    // Set coeff for implicit under-relaxation
    scalar relaxCoeff = 0.0;
    if (mesh_.relaxEquation(schemesField_))
    {
        relaxCoeff = mesh_.equationRelaxationFactor(schemesField_);
    }

    word divScheme("div(phi," + schemesField_ + ")");

    bool converged = false;

    volScalarField alphaClip = pos0(alpha_ - alphaCutOff_);

    for (label i = 0; i < nCorr_ ; i++)
    {
        Info << "Iteration: " << i << endl;

        if (normalize_)
        {
            normalizeA2();
        }

        // Compute eigenvectors
        closureModel_->updateEigenValsAndVecs();
        
        closureModel_->computeClosure(D_doubleDot_A4_, D_);

        volSymmTensorField A_HD_MRD = symm( 
                                        (W_ & A2_) - (A2_ & W_) 
                                        +xi_*(
                                                (D_ & A2_) + (A2_ & D_) - (2.0*D_doubleDot_A4_)
                                            )
                                    );

        computeC();

        closureModel_->computeClosure(D_doubleDot_A4_, C_);

        // // Original formulation
        // A_HD_MRD += 2.0*shrRate_*(C_ - tr(C_)*A2_);

        // Formulation by Favaloro and Tucker
        A_HD_MRD += shrRate_*(
                                2.0*C_ 
                              - 2.0*tr(C_)*A2_
                              - 5.0*symm((C_ & A2_) + (A2_ & C_))
                              + 10.0*D_doubleDot_A4_*dimensionedScalar("unitConversion", dimTime, 1.0)
                            );

        fvSymmTensorMatrix dA2dtEqn
        (
            fvm::ddt(A2_) + fvm::div(phi_, A2_, divScheme)
            ==
            alphaClip*(A_HD_MRD)
        );

        dA2dtEqn.relax(relaxCoeff);
        
        const scalar maxResidual = cmptMax(
                                            dA2dtEqn.solve(mesh_.solverDict(schemesField_)).initialResidual()
                                          );
        
        if (maxResidual < absTol_ && nCorr_ != 1)
        {
            Info << "Converged in: " << i + 1 <<" iterations" << endl;
            converged = true;
            i = nCorr_;
        }

        // Let eigen computation happen in the next iteration if needed
        closureModel_->updatedEigen(false);
    }

    if (!converged && nCorr_ != 1)
    {
        Info << "Did not converged in: " << nCorr_ <<" iterations" << endl;
    }
}

Foam::Ostream& Foam::fiberOrientation::operator<<
(
    Ostream& os,
    const MRD& ft
)
{
    ft.write(os);
    return os;
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fiberOrientation::MRD::computeC()
{
    const volTensorField& eigVecs = closureModel_->eigVecs();
    
    symmTensorField& C = C_.primitiveFieldRef();

    forAll(C, cellI)
    {
        // Compute e @ D @ e^T. Where D is the diagonal tensor with D1, D2 and D3 as diagonal elements
        const tensor& e = eigVecs[cellI];

        const scalar tmp0 = D1_*e.zx();
        const scalar tmp1 = D2_*e.yx();
        const scalar tmp2 = D3_*e.xx();
        const scalar tmp3 = tmp0*e.zy() + tmp1*e.yy() + tmp2*e.xy();
        const scalar tmp4 = tmp0*e.zz() + tmp1*e.yz() + tmp2*e.xz();
        const scalar tmp5 = D1_*e.zy()*e.zz() + D2_*e.yy()*e.yz() + D3_*e.xy()*e.xz();

        C[cellI].xx() = D1_*(e.zx()*e.zx()) + D2_*(e.yx()*e.yx()) + D3_*(e.xx()*e.xx());
        C[cellI].xy() = tmp3;
        C[cellI].xz() = tmp4;
        C[cellI].yy() = D1_*(e.zy()*e.zy()) + D2_*(e.yy()*e.yy()) + D3_*(e.xy()*e.xy());
        C[cellI].yz() = tmp5;
        C[cellI].zz() = D1_*(e.zz()*e.zz()) + D2_*(e.yz()*e.yz()) + D3_*(e.xz()*e.xz());

        C[cellI] *= CI_;
    }
}

bool Foam::fiberOrientation::MRD::read()
{
    fiberOrientationModel::read();
    return true;
}


void Foam::fiberOrientation::MRD::write(Ostream& os) const
{
    dictionary dict(typeName + "Properties");
    dict.add("CI", CI_);
    dict.add("D1", D1_);
    dict.add("D2", D2_);
    dict.add("D3", D3_);

    fiberOrientationModel::writeWithSubDict(os, dict);
}

// ************************************************************************* //
