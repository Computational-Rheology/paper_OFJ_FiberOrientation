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

#include "folgarTucker.H"
#include "addToRunTimeSelectionTable.H"
        
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{
    defineTypeNameAndDebug(folgarTucker, 0);

    addToRunTimeSelectionTable                                             
    (                                                                      
        fiberOrientationModel,      // Mother Class                                        
        folgarTucker,               // Current Class                                     
        dictionary                 // dictionary                                           
    );  
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fiberOrientation::folgarTucker::folgarTucker
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
    )
{
    Info << *this << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fiberOrientation::folgarTucker::solve()
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

        closureModel_->computeClosure(D_doubleDot_A4_, D_);

        fvSymmTensorMatrix dA2dtEqn
        (
            fvm::ddt(A2_) + fvm::div(phi_, A2_, divScheme)
            ==
            alphaClip*(
                        symm( 
                                (W_ & A2_) - (A2_ & W_) 
                                +xi_*(
                                        (D_ & A2_) + (A2_ & D_) - (2.0*D_doubleDot_A4_)
                                     )
                            )
                        + 2.0*CI_*shrRate_*(symmTensor::I - 3.0*A2_)
                    )
        );

        dA2dtEqn.relax(relaxCoeff);
        
        const scalar maxResidual = cmptMax(
                                            dA2dtEqn.solve(
                                                            mesh_.solverDict(schemesField_)
                                                          ).initialResidual()
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
    const folgarTucker& ft
)
{
    ft.write(os);
    return os;
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
bool Foam::fiberOrientation::folgarTucker::read()
{
    fiberOrientationModel::read();
    return true;
}


void Foam::fiberOrientation::folgarTucker::write(Ostream& os) const
{
    dictionary dict(typeName + "Properties");
    dict.add("CI", CI_);

    fiberOrientationModel::writeWithSubDict(os, dict);
}

// ************************************************************************* //
