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
    Foam::fiberOrientation::FT_RSC

Description
    Folgar-Tucker (FT) model with Reduced Strain Closure (RSC)

SourceFiles
    FT_RSC.C

Literature
    J. Wang, J. F. O’Gara, and C. L. Tucker, 
    “An objective model for slow orientation kinetics in concentrated fiber suspensions: Theory and rheological evidence,” 
    J. Rheol. (N. Y. N. Y)., vol. 52, no. 5, pp. 1179–1200, 2008, doi: 10.1122/1.2946437.


    J. H. Phelps and C. L. Tucker, 
    “An anisotropic rotary diffusion model for fiber orientation in short- and long-fiber thermoplastics,”
     J. Nonnewton. Fluid Mech., vol. 156, no. 3, pp. 165–176, 2009, doi: 10.1016/j.jnnfm.2008.08.002.

\*---------------------------------------------------------------------------*/

#ifndef fiberOrientationModel_FT_RSC_H
#define fiberOrientationModel_FT_RSC_H

#include "fiberOrientationModel.H"
#include "closureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{
    
class FT_RSC;

Ostream& operator<< (Ostream&, const FT_RSC&);

/*---------------------------------------------------------------------------*\
                           Class FT_RSC Declaration
\*---------------------------------------------------------------------------*/

class FT_RSC
:
    public fiberOrientationModel 
{
private:

    // Dictionary
    const dictionary& modelProperties_;

    // Scalars
    const scalar CI_;
    const scalar k_;

    // Functions 
    virtual bool read(); 
    virtual void write(Ostream& os) const; 


public:

    //- Runtime type information
    TypeName("FT-RSC");


    // Constructors

        //- Construct from components
        FT_RSC
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const volVectorField& U,
            const surfaceScalarField& phi         
        );

    //- Destructor
    virtual ~FT_RSC() = default;

    // Member Functions
        virtual void solve();

    // Ostream Operator
    friend Ostream& operator<< (Ostream& os, const FT_RSC& ft_rsc);   
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fiberOrientation
} // End namespace Foam



#endif

// ************************************************************************* //
