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
    Foam::fiberOrientation::closureModel

Description
    Abstract base class for closure models in fiber orientation modeling. 

SourceFiles
    closureModel.C

\*---------------------------------------------------------------------------*/

#ifndef closureModel_H
#define closureModel_H

#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{

/*---------------------------------------------------------------------------*\
                           Class closureModel Declaration
\*---------------------------------------------------------------------------*/

class closureModel
{
protected:

    // Protected data   
        // Fiber orientation tensor
        const volSymmTensorField& A2_;

        // boolean to prevent multiple computation of eigen quantities in one iteration
        bool updatedEig_;

        // eigenValues
        autoPtr<volVectorField> eigVals_;

        // eigenVectors
        autoPtr<volTensorField> eigVecs_;

        // Member Functions
        void createEigenValsAndVecs();
        void computeEigenValsAndVecs();
    
    //- No copy constructor
    closureModel(const closureModel&) = delete;

    //- No copy assignment
    void operator=(const closureModel&) = delete;

public:

    //- Runtime type information
    TypeName("closureModel");

    // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            closureModel,
            dictionary,
            (
                const dictionary& dict,
                const volSymmTensorField& A2
            ),
            (dict, A2)
        );

    // Selectors
        //- Return a reference to the selected fiber model
        static autoPtr<closureModel> New
        (
            const dictionary& dict,
            const volSymmTensorField& A2
        );

    // Constructors

        //- Construct from components
        closureModel
        (
            const dictionary& dict,
            const volSymmTensorField& A2
        );

    //- Destructor
    virtual ~closureModel() = default;

    // Member Functions
        
    // Virtual functions
        virtual void computeClosure(
                                        volSymmTensorField& D_doubleDot_A4,
                                        const volSymmTensorField& D
                                    ) = 0;

        virtual void computeRSCClosure(
                                        volSymmTensorField& D_doubleDot_A4, 
                                        const volSymmTensorField& D, 
                                        const scalar& k
                                        );


    // Member Functions
        void updateEigenValsAndVecs(bool forceCalculation=true);

        // Access
        const volVectorField& eigVals() const
        {
            return eigVals_();
        }

        const volTensorField& eigVecs() const
        {
            return eigVecs_();
        }

        void updatedEigen(bool flag)
        {
            updatedEig_ = flag;
        }
};

} // End namespace fiberOrientation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
