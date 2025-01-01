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
    Class to manage retarding principal rate (RPR). 
    Used in iARD-RPR and pARD-RPR

\*---------------------------------------------------------------------------*/

#ifndef fiberOrientationModel_RPR_H
#define fiberOrientationModel_RPR_H

#include "closureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{


/*---------------------------------------------------------------------------*\
                           Class RPR Declaration
\*---------------------------------------------------------------------------*/

class RPR
{
protected:
    // Closure model
    const closureModel& closureModel_;

    // Scalars
    const scalar k_;
    const scalar alpha_;
    const scalar beta_;

    // Fields
    volSymmTensorField Adot_RPR_;

public:

    //- Runtime type information
    TypeName("RPR");

    // Constructors

        //- Construct from components
        RPR
        (
            const closureModel& closureModel,
            const dictionary& modelProperties,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~RPR() = default;

    // Functions
    void computeRPR(const volSymmTensorField& Adot_HD_iARD);
    void computeLambdaDotIOK(diagTensor& store, const diagTensor& lambda_Adot_HR_iARD);

    scalar k() const {return k_;}
    scalar alpha() const {return alpha_;}
    scalar beta() const {return beta_;}
    const volSymmTensorField& Adot_RPR() const {return Adot_RPR_;}
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fiberOrientation
} // End namespace Foam



#endif

// ************************************************************************* //
