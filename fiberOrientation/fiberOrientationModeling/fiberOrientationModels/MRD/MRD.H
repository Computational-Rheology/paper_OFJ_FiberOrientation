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
    Foam::fiberOrientation::MRD

Description
    Moldflow Rotary diffusion model

SourceFiles
    MRD.C
    
Literature
    A. Bakharev, H. Yu, S. Ray, R. Speight, and J. Wang, 
    “Using new anisotropic rotational diffusion model to improve prediction of short fibers in thermoplastic injection molding,”
    Annu. Tech. Conf. - ANTEC, Conf. Proc., vol. 2018-May, 2018.

    A. J. Favaloro and C. L. Tucker, 
    “Analysis of anisotropic rotary diffusion models for fiber orientation,”
    Compos. Part A Appl. Sci. Manuf., vol. 126, no. August, 2019, doi: 10.1016/j.compositesa.2019.105605. 
    *Note: Uses a different formulation for this model (Phelps-Tucker formulation)

    S. K. Kugler, G. M. Lambert, C. Cruz, A. Kech, T. A. Osswald, and D. G. Baird,
    “Macroscopic fiber orientation model evaluation for concentrated short fiber reinforced polymers in comparison to experimental data,” 
    Polym. Compos., vol. 41, no. 7, pp. 2542–2556, 2020, doi: 10.1002/pc.25553.

\*---------------------------------------------------------------------------*/

#ifndef fiberOrientationModel_MRD_H
#define fiberOrientationModel_MRD_H

#include "fiberOrientationModel.H"
#include "closureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{

class MRD;

Ostream& operator<< (Ostream&, const MRD&);

/*---------------------------------------------------------------------------*\
                           Class MRD Declaration
\*---------------------------------------------------------------------------*/

class MRD
:
    public fiberOrientationModel 
{
private:

    // Dictionary
    const dictionary& modelProperties_;
    
    // Scalars
    const scalar CI_;
    const scalar D1_;
    const scalar D2_;
    const scalar D3_;

    volSymmTensorField C_;

    // Member Functions 
    virtual void write(Ostream& os) const;
    virtual bool read();

    // Compute anisotropic rotary diffusion
    void computeC();

public:

    //- Runtime type information
    TypeName("MRD");


    // Constructors

        //- Construct from components
        MRD
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const volVectorField& U,
            const surfaceScalarField& phi
        );

    //- Destructor
    virtual ~MRD() = default;
    
    // Member Functions
        virtual void solve();

    // Ostream Operator
        friend Ostream& operator<< (Ostream& os, const MRD& ft);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fiberOrientation
} // End namespace Foam



#endif

// ************************************************************************* //
