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

#include "iARD_RPR.H"
#include "addToRunTimeSelectionTable.H"
        
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{
    defineTypeNameAndDebug(iARD_RPR, 0);

    addToRunTimeSelectionTable                                             
    (                                                                      
        fiberOrientationModel,      // Mother Class                                        
        iARD_RPR,                   // Current Class                                     
        dictionary                  // dictionary                                           
    );  
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fiberOrientation::iARD_RPR::iARD_RPR
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
CM_(readScalar(modelProperties_.lookup("CM"))),
Dsquare_
(
    IOobject
    (
        "fiberModel_Dsquare",
        U.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimensionedSymmTensor(dimless/(dimTime * dimTime), Foam::Zero)
),
Dr_
(
    IOobject
    (
        "fiberModel_Dr",
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
    Dsquare_ = symm(D_ & D_);

    const symmTensorField& Dsquare = Dsquare_.primitiveField();

    Dr_.primitiveFieldRef() = CI_*(symmTensor::I - CM_*(Dsquare/Foam::sqrt(0.5*(Dsquare && Dsquare))));

    Info << *this << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fiberOrientation::iARD_RPR::solve()
{   
    correctFlowFieldTensors();

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

        volSymmTensorField Adot_HD_iARD = symm(  
                                                ((W_ & A2_) - (A2_ & W_))
                                                + xi_*((D_ & A2_) + (A2_ & D_)-(2.0*D_doubleDot_A4_))
                                              );

        closureModel_->computeClosure(D_doubleDot_A4_, Dr_);

        Adot_HD_iARD += shrRate_*(
                                    2.0*Dr_ 
                                    - 2.0*tr(Dr_)*A2_
                                    - 5.0*symm((Dr_ & A2_) + (A2_ & Dr_))
                                    + 10.0*D_doubleDot_A4_*dimensionedScalar("unitConversion", dimTime, 1.0)
                                );

        // Compute RPR
        RPR_.computeRPR(Adot_HD_iARD);

        fvSymmTensorMatrix dA2dtEqn
        (
            fvm::ddt(A2_) + fvm::div(phi_, A2_, divScheme)
            ==
            alphaClip*(Adot_HD_iARD + RPR_.Adot_RPR())
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
    const iARD_RPR& iard_rpr
)
{
    iard_rpr.write(os);
    return os;
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fiberOrientation::iARD_RPR::read()
{
    fiberOrientationModel::read();
    return true;
}


void Foam::fiberOrientation::iARD_RPR::correctFlowFieldTensors()
{
    if (updateFlowFields_)
    {
        fiberOrientationModel::correctFlowFieldTensors();
        
        Dsquare_ = symm(D_ & D_);

        Dr_ = CI_*(symmTensor::I - CM_*(Dsquare_/Foam::sqrt(0.5*(Dsquare_ && Dsquare_))));
    }
}

void Foam::fiberOrientation::iARD_RPR::write(Ostream& os) const
{
   dictionary dict(typeName + "Properties");
   dict.add("CI", CI_);
   dict.add("CM", CM_);
   dict.add("alpha", RPR_.alpha());
   dict.add("beta", RPR_.beta());

    fiberOrientationModel::writeWithSubDict(os, dict);
}

// ************************************************************************* //
