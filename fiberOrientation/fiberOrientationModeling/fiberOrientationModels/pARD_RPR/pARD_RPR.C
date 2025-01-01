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

#include "pARD_RPR.H"
#include "addToRunTimeSelectionTable.H"
        
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{
    defineTypeNameAndDebug(pARD_RPR, 0);

    addToRunTimeSelectionTable                                             
    (                                                                      
        fiberOrientationModel,      // Mother Class                                        
        pARD_RPR,                   // Current Class                                     
        dictionary                 // dictionary                                           
    );  
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fiberOrientation::pARD_RPR::pARD_RPR
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
        ? readScalar(modelProperties_.lookup("CI")) 
        : computeCI() 
    ),
omega_(modelProperties_.get<scalar>("omega")),
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
),
RPR_(this->closureModel_(), modelProperties_, mesh_)
{
    Info << *this << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fiberOrientation::pARD_RPR::solve()
{   
    fiberOrientationModel::correctFlowFieldTensors();

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

        // Compute eigen-quantities                                           
        closureModel_->updateEigenValsAndVecs();

        closureModel_->computeClosure(D_doubleDot_A4_, D_);

        volSymmTensorField Adot_HD_pARD = symm(  
                                                ((W_ & A2_) - (A2_ & W_))
                                                + xi_*((D_ & A2_) + (A2_ & D_)-(2.0*D_doubleDot_A4_))
                                              );

        computeC();

        closureModel_->computeClosure(D_doubleDot_A4_, C_);

        Adot_HD_pARD += shrRate_*(
                                    2.0*C_ 
                                    - 2.0*tr(C_)*A2_
                                    - 5.0*symm((C_ & A2_) + (A2_ & C_))
                                    + 10.0*D_doubleDot_A4_*dimensionedScalar("unitConversion", dimTime, 1.0)
                                );
                                
        RPR_.computeRPR(Adot_HD_pARD);

        fvSymmTensorMatrix dA2dtEqn
        (
            fvm::ddt(A2_) + fvm::div(phi_, A2_, divScheme)
            ==
            alphaClip*(Adot_HD_pARD + RPR_.Adot_RPR())
        );

        dA2dtEqn.relax(relaxCoeff);
        
        scalar maxResidual = cmptMax(
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
    const pARD_RPR& ft
)
{
    ft.write(os);
    return os;
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fiberOrientation::pARD_RPR::computeC()
{
    const volTensorField& eigVecs = closureModel_->eigVecs();
    
    symmTensorField& C = C_.primitiveFieldRef();

    forAll(C, cellI)
    {
        // Compute e @ D @ e^T. Where D is the diagonal tensor with D1,D2 and D3 as diagonal elements
        const tensor& e = eigVecs[cellI];

        const scalar tmp0 = 1 - omega_;
        const scalar tmp1 = omega_*e.yx();
        const scalar tmp2 = tmp0*e.xx();
        const scalar tmp3 = tmp1*e.yy() + tmp2*e.xy() + e.zx()*e.zy();
        const scalar tmp4 = tmp1*e.yz() + tmp2*e.xz() + e.zx()*e.zz();
        const scalar tmp5 = omega_*e.yy()*e.yz() + tmp0*e.xy()*e.xz() + e.zy()*e.zz();

        C[cellI].xx() = omega_*(e.yx()*e.yx()) + tmp0*(e.xx()*e.xx()) + e.zx()*e.zx();
        C[cellI].xy() = tmp3;
        C[cellI].xz() = tmp4;
        C[cellI].yy() = omega_*(e.yy()*e.yy()) + tmp0*(e.xy()*e.xy()) + e.zy()*e.zy();
        C[cellI].yz() = tmp5;
        C[cellI].zz() = omega_*(e.yz()*e.yz()) + tmp0*(e.xz()*e.xz()) + e.zz()*e.zz();

        C[cellI] *= CI_;
    }
}

bool Foam::fiberOrientation::pARD_RPR::read()
{
    fiberOrientationModel::read();
    return true;
}


void Foam::fiberOrientation::pARD_RPR::write(Ostream& os) const
{
    dictionary dict(typeName + "Properties");
    dict.add("CI", CI_);
    dict.add("omega", omega_);
    dict.add("alpha", RPR_.alpha());
    dict.add("beta", RPR_.beta());

    fiberOrientationModel::writeWithSubDict(os, dict);
}

// ************************************************************************* //
