/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "fiberOrientationModeling.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fiberOrientationModeling, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        fiberOrientationModeling,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fiberOrientationModeling::fiberOrientationModeling
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phaseName_(dict.getOrDefault<word>("phase", "none")),
    phiName_(   phaseName_ == "none" 
                ? 
                dict.getOrDefault<word>("phi", "phi"):
                dict.getOrDefault<word>("phasePhiCompressed", "alphaPhiUn") ),
    U_( mesh_.lookupObject<volVectorField>("U") ),
    phi_( mesh_.lookupObject<surfaceScalarField>(phiName_) ),
    fiberOrientationModel_(fiberOrientation::fiberOrientationModel::New(dict, mesh_, U_, phi_))
{
    read(dict);

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];
        if (isA<emptyPolyPatch>(pp))
        {
            FatalErrorInFunction
                << "Fiber Orientation modeling is not ready for 2D cases!" << exit(FatalError);
        }
    }

    Info << nl
         << "FO is created" << nl
         << "Phase name is:\t" << phaseName_ << nl
         << "Phi name is:\t" << phiName_ << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fiberOrientationModeling::~fiberOrientationModeling()
{
    Info << "FiberOrientation FO was destroyed" << endl;
} 


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fiberOrientationModeling::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("phase", phaseName_);

    if (phaseName_ != "none")
    {
        phiName_ = dict.getOrDefault<word>("phasePhiCompressed", "alphaPhiUn");
    }
    else
    {
        phiName_ = dict.getOrDefault<word>("phi", "phi");
    }

    return true;
}


bool Foam::functionObjects::fiberOrientationModeling::execute()
{
    Log << type() << ":\tSolving fiber orientation:" << endl;

    fiberOrientationModel_->solve();

    Log << endl;
    
    return true;
}

bool Foam::functionObjects::fiberOrientationModeling::write()
{
    return true;
}


// ************************************************************************* //
