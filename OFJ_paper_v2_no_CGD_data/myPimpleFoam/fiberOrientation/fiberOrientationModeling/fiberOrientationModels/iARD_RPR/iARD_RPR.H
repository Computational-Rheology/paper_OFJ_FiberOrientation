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
    Foam::fiberOrientation::iARD_RPR

Description
    iARD-RPR model

SourceFiles
    iARD_RPR.C

Literature
    H.-C. Tseng, R.-Y. Chang, and C.-H. Hsu, 
    “An objective tensor to predict anisotropic fiber orientation in concentrated suspensions,” 
    J. Rheol. (N. Y. N. Y)., vol. 60, no. 2, pp. 215–224, 2016, doi: https://doi.org/10.1122/1.4939098.

\*---------------------------------------------------------------------------*/

#ifndef fiberOrientationModel_iARD_RPR_H
#define fiberOrientationModel_iARD_RPR_H

#include "fiberOrientationModel.H"
#include "closureModel.H"
#include "RPR.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{

class iARD_RPR;

Ostream& operator<< (Ostream&, const iARD_RPR&);

/*---------------------------------------------------------------------------*\
                           Class iARD_RPR Declaration
\*---------------------------------------------------------------------------*/

class iARD_RPR
:
    public fiberOrientationModel
{
    private:
        // Dictionary
        const dictionary& modelProperties_;
        
        // Scalars
        const scalar CI_;
        const scalar CM_;

        // Fields
            volSymmTensorField Dsquare_; 
            volSymmTensorField Dr_; 

        // RPR 
            RPR RPR_;

        // Functions
            virtual bool read();
            virtual void write(Ostream& os) const;
            virtual void correctFlowFieldTensors();

public:

    //- Runtime type information
    TypeName("iARD_RPR");

    // Constructors

        //- Construct from components
        iARD_RPR
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const volVectorField& U,
            const surfaceScalarField& phi
        );

    //- Destructor
    virtual ~iARD_RPR() = default;
    
    // Member Functions
    
        virtual void solve();
        
        // Ostream Operator
        friend Ostream& operator<< (Ostream& os,  const iARD_RPR& iard_rpr);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fiberOrientation
} // End namespace Foam

#endif

// ************************************************************************* //
