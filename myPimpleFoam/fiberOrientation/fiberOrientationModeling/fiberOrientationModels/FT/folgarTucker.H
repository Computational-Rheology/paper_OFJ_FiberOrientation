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
    Foam::fiberOrientation::folgarTucker

Description
    Folgar-Tucker (FT) model

SourceFiles
    folgarTucker.C

Literature
    S. G. Advani and C. L. Tucker, 
    “The Use of Tensors to Describe and Predict Fiber Orientation in Short Fiber Composites,” 
    J. Rheol. (N. Y. N. Y)., vol. 31, no. 8, pp. 751–784, 1987, doi: 10.1122/1.549945.

\*---------------------------------------------------------------------------*/

#ifndef fiberOrientationModel_folgarTucker_H
#define fiberOrientationModel_folgarTucker_H

#include "fiberOrientationModel.H"
#include "closureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{

class folgarTucker;

Ostream& operator<< (Ostream&, const folgarTucker&);

/*---------------------------------------------------------------------------*\
                           Class folgarTucker Declaration
\*---------------------------------------------------------------------------*/

class folgarTucker
:
    public fiberOrientationModel 
{
private:

    // Dictionary
    const dictionary& modelProperties_;
    
    // Scalars
    const scalar CI_;

    // Functions 
    virtual bool read();
    virtual void write(Ostream& os) const;

public:

    //- Runtime type information
    TypeName("folgarTucker");


    // Constructors

        //- Construct from components
        folgarTucker
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const volVectorField& U,
            const surfaceScalarField& phi
        );

    //- Destructor
    virtual ~folgarTucker() = default;
    
    // Member Functions
        virtual void solve();

    // Ostream Operator
        friend Ostream& operator<< (Ostream& os, const folgarTucker& ft);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fiberOrientation
} // End namespace Foam



#endif

// ************************************************************************* //
