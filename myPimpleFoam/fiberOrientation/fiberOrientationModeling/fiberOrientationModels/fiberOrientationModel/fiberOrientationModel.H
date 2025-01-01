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
    Foam::fiberOrientationModel

Description
    Abstract base class for orientation models. 

SourceFiles
    fiberOrientationModel.C

\*---------------------------------------------------------------------------*/

#ifndef fiberOrientationModel_H
#define fiberOrientationModel_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "closureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation 
{

/*---------------------------------------------------------------------------*\
                           Class fiberOrientationModel Declaration
\*---------------------------------------------------------------------------*/

class fiberOrientationModel;

Ostream& operator<<
(
    Ostream&,
    const fiberOrientationModel&
);


class fiberOrientationModel
{
private:

    scalar computeXi() const;
    scalar fiberVolumeFraction() const;
    volSymmTensorField& createA2(); 
    const volScalarField& constLookupOrConstruct(const word& name);  

protected:

    // Protected data
        const fvMesh&               mesh_;

        // Dictionary
        dictionary                  fiberPropertiesDict_;

        // Fields
        const volVectorField&       U_;       
        const surfaceScalarField&   phi_; 
        const volScalarField&       alpha_;
        
        volSymmTensorField&         A2_; 
        volTensorField              L_;      
        volSymmTensorField          D_;  
        volTensorField              W_;
        volScalarField              shrRate_;           
        volSymmTensorField          D_doubleDot_A4_;  
      
        // Scalars
        const scalar xi_;
        const scalar alphaCutOff_; // (optional)

        // Closure model
        autoPtr<fiberOrientation::closureModel> closureModel_;

        // Functions
        scalar computeCI() const;
        const dictionary& fiberDict() const;
        void writeWithSubDict(Ostream& os, const dictionary& subDict) const;
            
            // Virtual functions
            virtual void normalizeA2();              
            virtual void correctFlowFieldTensors();
            
            // Pure virtual functions
            virtual bool read() = 0; 
            virtual void write(Ostream& os) const;  
            

        // Solution Properties
        dictionary     solutionPropertiesDict_;

        //- Number of corrector iterations (optional)
        label nCorr_;

        //- Name of field whose schemes are used (optional)
        word schemesField_;   

        // Absolute residual    (optional)
        scalar absTol_;

        // Choose to normalize A2   (optional)
        Switch normalize_;

        // Choose to update fields  (optional)
        Switch updateFlowFields_;

    //- No copy constructor
    fiberOrientationModel(const fiberOrientationModel&) = delete;

    //- No copy assignment
    void operator=(const fiberOrientationModel&) = delete;

public:

    //- Runtime type information
    TypeName("fiberOrientationModel");

    // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            fiberOrientationModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (dict, mesh, U, phi)
        );

    // Selectors
        //- Return a reference to the selected fiber model
        static autoPtr<fiberOrientationModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const volVectorField& U,
            const surfaceScalarField& phi
        );

    // Constructors

        //- Construct from components
        fiberOrientationModel
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const volVectorField& U,
            const surfaceScalarField& phi
        );

    //- Destructor
    virtual ~fiberOrientationModel() = default;

    // Member Functions
        
        // Pure virtual functions
        virtual void solve() = 0;

        // Ostream Operator
        friend Ostream& operator<<(Ostream& os, const fiberOrientationModel& fm);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace fiberOrientation
} // End namespace Foam


#endif

// ************************************************************************* //
