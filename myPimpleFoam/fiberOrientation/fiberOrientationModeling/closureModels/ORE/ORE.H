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


Description
    The following function was generated with symbolic computation

Literature
    B. E. VerWeyst.
    "Numerical predictions of flow-induced fiber orientation in 3-D geometries". PhD thesis,
    University of Illinois at Urbana-Champaign, Urbana, IL, 1998.
    
\*-------------------------------------------------------------------------------------------------------*/

#ifndef closureModel_OREClosure_H
#define closureModel_OREClosure_H

#include "closureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{
namespace closureModels
{
/*---------------------------------------------------------------------------*\
                           Class ORE Declaration
\*---------------------------------------------------------------------------*/

class ORE
:
    public closureModel
{
protected:

public:

    //- Runtime type information
    TypeName("ORE");

    // Constructors

        //- Construct from components
        ORE
        (
            const dictionary& dict,
            const volSymmTensorField& A2
        );

    //- Destructor
    virtual ~ORE() = default;

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

        void D_doubleDot_A4_ORE
        (
            symmTensor& result, 
            const symmTensor& D,
            const Foam::vector& eigenValue,
            const Foam::tensor& eigenVector
        );

        void RSC_D_doubleDot_A4_ORE
        (
            symmTensor& result,
            const symmTensor& D,
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
