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

#include "IBOF.H"
#include "addToRunTimeSelectionTable.H"
        
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{   
namespace fiberOrientation
{
namespace closureModels
{
    defineTypeNameAndDebug(IBOF, 0);

    addToRunTimeSelectionTable
    (
        closureModel, 
        IBOF, 
        dictionary
    );  
}
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fiberOrientation::closureModels::IBOF::IBOF
( 
    const dictionary& dict,
    const volSymmTensorField& A2
)
:
closureModel(dict, A2)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "D_doubleDot_A4_IBOF.H"
#include "RSC_D_doubleDot_A4_IBOF.H"

void Foam::fiberOrientation::closureModels::IBOF::computeClosure
(
    volSymmTensorField& D_doubleDot_A4,
    const volSymmTensorField& D
)
{
    forAll(D_doubleDot_A4, cellI)
    {      
        D_doubleDot_A4_IBOF( D_doubleDot_A4[cellI], D[cellI], A2_[cellI]);
    }
}


void Foam::fiberOrientation::closureModels::IBOF::computeRSCClosure
(
    volSymmTensorField& D_doubleDot_A4, 
    const volSymmTensorField& D, 
    const scalar& k
)
{
    updateEigenValsAndVecs(false);

    const volVectorField& eigVals = eigVals_.ref();
    const volTensorField& eigVecs = eigVecs_.ref();

    forAll(D_doubleDot_A4, cellI)
    {   
        RSC_D_doubleDot_A4_IBOF(
                                    D_doubleDot_A4[cellI], 
                                    D[cellI], 
                                    A2_[cellI],
                                    eigVals[cellI],
                                    eigVecs[cellI], 
                                    k
                                );
    }
}
        
// ************************************************************************* //
