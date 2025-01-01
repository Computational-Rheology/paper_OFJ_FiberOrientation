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
    Foam::fiberOrientation::pARD_RPR

Description
    principal Anisotropic Rotary Diffusion model with RPR

SourceFiles
    pARD_RPR.C


Literature
    H.-C. Tseng, R.-Y. Chang, and C.-H. Hsu, 
    “The use of principal spatial tensor to predict anisotropic fiber orientation in concentrated fiber suspensions,”
    J. Rheol. (N. Y. N. Y)., vol. 62, no. 1, pp. 313–320, 2018, doi: 10.1122/1.4998520.

\*---------------------------------------------------------------------------*/

#ifndef fiberOrientationModel_pARD_RPR_H
#define fiberOrientationModel_pARD_RPR_H

#include "iARD_RPR.H"
#include "closureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{

class pARD_RPR;

Ostream& operator<< (Ostream&, const pARD_RPR&);

/*---------------------------------------------------------------------------*\
                           Class pARD_RPR Declaration
\*---------------------------------------------------------------------------*/

class pARD_RPR
:
    public fiberOrientationModel 
{
private:

    // Model dictionary
    const dictionary& modelProperties_;

    // scalars
    const scalar CI_;
    const scalar omega_;
    volSymmTensorField C_;

    // RPR component
    RPR RPR_;

    // Member Functions 
    virtual void write(Ostream& os) const;
    virtual bool read();

    // Compute anisotropic rotary diffusion
    void computeC();

public:

    //- Runtime type information
    TypeName("pARD_RPR");


    // Constructors

        //- Construct from components
        pARD_RPR
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const volVectorField& U,
            const surfaceScalarField& phi
        );

    //- Destructor
    virtual ~pARD_RPR() = default;
    
    // Member Functions
        virtual void solve();

    // Ostream Operator
        friend Ostream& operator<< (Ostream& os, const pARD_RPR& ft);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fiberOrientation
} // End namespace Foam



#endif

// ************************************************************************* //
