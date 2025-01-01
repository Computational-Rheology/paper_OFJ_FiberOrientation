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

#include "closureModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{
    defineTypeNameAndDebug(closureModel, 0);

    defineRunTimeSelectionTable 
    ( 
        closureModel,
        dictionary
    ); 
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::fiberOrientation::closureModel::closureModel
(
    const dictionary& dict,
    const volSymmTensorField& A2
)
:
A2_(A2),
updatedEig_(false),
eigVals_(nullptr),
eigVecs_(nullptr)
{}


void Foam::fiberOrientation::closureModel::computeRSCClosure
(
    volSymmTensorField& D_doubleDot_A4,
    const volSymmTensorField& D,
    const scalar& k
)
{
    NotImplemented;
}


void Foam::fiberOrientation::closureModel::createEigenValsAndVecs()
{
    if(!eigVals_)
    {
        eigVals_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "eigenValues",
                    A2_.time().timeName(),
                    A2_.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                eigenValues(A2_)
            )
        );
    }

    if(!eigVecs_)
    {
        eigVecs_.reset
        (
            new volTensorField
            (
                IOobject
                (
                    "eigenVectors",
                    A2_.time().timeName(),
                    A2_.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false 
                ),
                eigenVectors(A2_)
            )
        );
    }
}

void Foam::fiberOrientation::closureModel::computeEigenValsAndVecs()
{
    if(!eigVals_ || !eigVecs_)
    {
        FatalErrorInFunction
            << "No eigenValues or eigenVectors are created!" << exit(FatalError);
    }

    volVectorField& eigVals = eigVals_.ref();
    volTensorField& eigVecs = eigVecs_.ref();

    forAll(eigVals, cellI)
    {
        eigVals[cellI] =  eigenValues(A2_[cellI]);
        eigVecs[cellI]  = eigenVectors(A2_[cellI]);
    }
}


void Foam::fiberOrientation::closureModel::updateEigenValsAndVecs
(
    bool forceCalculation
)
{
    if(!eigVals_ || !eigVecs_)
    {
        createEigenValsAndVecs();
        updatedEig_= true;
        return;
    }

    if (!updatedEig_ || forceCalculation)
    {
        computeEigenValsAndVecs();
        updatedEig_ = true;
    }
}

// ************************************************************************* //
 