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
    Foam::fiberOrientation::closureModels::IBOF

Description
    IBOF closure. See Literature.

SourceFiles
    IBOF.C

Literature
    D. H. Chung and T. H. Kwon, 
    “Invariant-based optimal fitting closure approximation for the numerical prediction of flow-induced fiber orientation,” 
    J. Rheol. (N. Y. N. Y)., vol. 46, no. 1, pp. 169–194, 2002, doi: 10.1122/1.1423312.


\*---------------------------------------------------------------------------*/

#ifndef closureModel_IBOFClosure_H
#define closureModel_IBOFClosure_H

#include "closureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{
namespace closureModels
{
/*---------------------------------------------------------------------------*\
                           Class IBOF Declaration
\*---------------------------------------------------------------------------*/

class IBOF
:
    public closureModel
{

public:

    //- Runtime type information
    TypeName("IBOF");

    // Constructors

        //- Construct from components
        IBOF
        (
            const dictionary& dict,
            const volSymmTensorField& A2
        );

    //- Destructor
    virtual ~IBOF() = default;

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

        void D_doubleDot_A4_IBOF
        (
            symmTensor& result, 
            const symmTensor& D, 
            const symmTensor& A
        );

        void RSC_D_doubleDot_A4_IBOF
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
