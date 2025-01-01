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

Class
    Foam::fiberOrientation::closureModels::hybrid

Description
    Hybrid closure. See Literature.

SourceFiles
    hybrid.C

Literature
    S. G. Advani and C. L. Tucker,
    “Closure approximations for three‐dimensional structure tensors,” 
    J. Rheol. (N. Y. N. Y)., vol. 34, no. 3, pp. 367–386, 1990, 
    doi: 10.1122/1.550133.
    
\*---------------------------------------------------------------------------*/

#ifndef closureModel_hybridClosure_H
#define closureModel_hybridClosure_H

#include "closureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{
namespace closureModels
{
/*---------------------------------------------------------------------------*\
                           Class hybrid Declaration
\*---------------------------------------------------------------------------*/

class hybrid
:
    public closureModel 
{

public:

    //- Runtime type information
    TypeName("hybrid");

    // Constructors

        //- Construct from components
        hybrid
        (
            const dictionary& dict,
            const volSymmTensorField& A2
        );

    //- Destructor
    virtual ~hybrid() = default;

    // Member Functions

        virtual void computeClosure(
                                        volSymmTensorField& D_doubleDot_A4,
                                        const volSymmTensorField& D
                                    );

        virtual void computeRSCClosure(
                                        volSymmTensorField& D_doubleDot_A4, 
                                        const volSymmTensorField& D, 
                                        const scalar& k
                                        );        

        void D_doubleDot_A4_hybrid
        (
            symmTensor& result, 
            const symmTensor& D, 
            const symmTensor& A
        );

        void RSC_D_doubleDot_A4_hybrid
        (
            symmTensor& result,
            const symmTensor& D,
            const symmTensor& A,
            const vector& eigenValue,
            const tensor& eigenVector,
            const scalar& k
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace closureModels
} // End namespace fiberOrientation
} // End namespace Foam



#endif

// ************************************************************************* //
