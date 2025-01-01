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

#include "RPR.H"
        
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fiberOrientation
{
    defineTypeNameAndDebug(RPR, 0);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fiberOrientation::RPR::RPR
( 
    const closureModel& closureModel,
    const dictionary& modelProperties,
    const fvMesh& mesh
)
:
closureModel_(closureModel),
k_(modelProperties.getOrDefault<scalar>("k", 0.0)),
alpha_(
            modelProperties.found("alpha") 
            ? readScalar(modelProperties.lookup("alpha")) 
            : 1.0 - readScalar(modelProperties.lookup("k"))
      ),
beta_(modelProperties.getOrDefault<scalar>("beta", 0.0)),
Adot_RPR_
(
    IOobject
    (
        "fiberModel_Adot_RPR",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimensionedSymmTensor(dimless/dimTime, Foam::Zero)
)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fiberOrientation::RPR::computeRPR(const volSymmTensorField& Adot_HD_iARD)
{
    const volTensorField& eigVecs = closureModel_.eigVecs();

    forAll(Adot_RPR_, cellI)
    {    
        const tensor& e = eigVecs[cellI];

        diagTensor lambda_Adot_HR_iARD;

        {
            const scalar tmp0 = Adot_HD_iARD[cellI].xx()*e.zx() + Adot_HD_iARD[cellI].xy()*e.zy() + Adot_HD_iARD[cellI].xz()*e.zz();
            const scalar tmp1 = Adot_HD_iARD[cellI].xy()*e.zx() + Adot_HD_iARD[cellI].yy()*e.zy() + Adot_HD_iARD[cellI].yz()*e.zz();
            const scalar tmp2 = Adot_HD_iARD[cellI].xz()*e.zx() + Adot_HD_iARD[cellI].yz()*e.zy() + Adot_HD_iARD[cellI].zz()*e.zz();
            const scalar tmp3 = Adot_HD_iARD[cellI].xx()*e.yx() + Adot_HD_iARD[cellI].xy()*e.yy() + Adot_HD_iARD[cellI].xz()*e.yz();
            const scalar tmp4 = Adot_HD_iARD[cellI].xy()*e.yx() + Adot_HD_iARD[cellI].yy()*e.yy() + Adot_HD_iARD[cellI].yz()*e.yz();
            const scalar tmp5 = Adot_HD_iARD[cellI].xz()*e.yx() + Adot_HD_iARD[cellI].yz()*e.yy() + Adot_HD_iARD[cellI].zz()*e.yz();
            const scalar tmp6 = Adot_HD_iARD[cellI].xx()*e.xx() + Adot_HD_iARD[cellI].xy()*e.xy() + Adot_HD_iARD[cellI].xz()*e.xz();
            const scalar tmp7 = Adot_HD_iARD[cellI].xy()*e.xx() + Adot_HD_iARD[cellI].yy()*e.xy() + Adot_HD_iARD[cellI].yz()*e.xz();
            const scalar tmp8 = Adot_HD_iARD[cellI].xz()*e.xx() + Adot_HD_iARD[cellI].yz()*e.xy() + Adot_HD_iARD[cellI].zz()*e.xz();

            lambda_Adot_HR_iARD.xx() = tmp0*e.zx() + tmp1*e.zy() + tmp2*e.zz();
            lambda_Adot_HR_iARD.yy() = tmp3*e.yx() + tmp4*e.yy() + tmp5*e.yz();
            lambda_Adot_HR_iARD.zz() = tmp6*e.xx() + tmp7*e.xy() + tmp8*e.xz();
        }

        diagTensor lambdaDot_IOK(Foam::Zero);

        computeLambdaDotIOK(lambdaDot_IOK, lambda_Adot_HR_iARD);

        {
            const scalar tmp0 = e.xx()*lambdaDot_IOK.zz();
            const scalar tmp1 = e.yx()*lambdaDot_IOK.yy();
            const scalar tmp2 = e.zx()*lambdaDot_IOK.xx();
            const scalar tmp3 = -tmp0*e.xy() - tmp1*e.yy() - tmp2*e.zy();
            const scalar tmp4 = -tmp0*e.xz() - tmp1*e.yz() - tmp2*e.zz();
            const scalar tmp5 = -e.xy()*e.xz()*lambdaDot_IOK.zz() - e.yy()*e.yz()*lambdaDot_IOK.yy() - e.zy()*e.zz()*lambdaDot_IOK.xx();

            Adot_RPR_[cellI].xx() = -lambdaDot_IOK.xx()*e.zx()*e.zx() - lambdaDot_IOK.yy()*e.yx()*e.yx() - lambdaDot_IOK.zz()*e.xx()*e.xx();
            Adot_RPR_[cellI].xy() = tmp3;
            Adot_RPR_[cellI].xz() = tmp4;
            Adot_RPR_[cellI].yy() = -lambdaDot_IOK.xx()*e.zy()*e.zy() - lambdaDot_IOK.yy()*e.yy()*e.yy() - lambdaDot_IOK.zz()*e.xy()*e.xy();
            Adot_RPR_[cellI].yz() = tmp5;
            Adot_RPR_[cellI].zz() = -lambdaDot_IOK.xx()*e.zz()*e.zz() - lambdaDot_IOK.yy()*e.yz()*e.yz() - lambdaDot_IOK.zz()*e.xz()*e.xz();
        }
    }
}

void Foam::fiberOrientation::RPR::computeLambdaDotIOK
(
    diagTensor& LambdaDot_IOK, 
    const diagTensor& lambda_Adot_HR_iARD
)
{
    LambdaDot_IOK.xx() = alpha_*(
                                lambda_Adot_HR_iARD[0]
                                -beta_*(lambda_Adot_HR_iARD[0]*lambda_Adot_HR_iARD[0]
                                        + 2.0*lambda_Adot_HR_iARD[1]*lambda_Adot_HR_iARD[2]) 
                            );

    LambdaDot_IOK.yy() = alpha_*(
                                lambda_Adot_HR_iARD[1]
                                    -beta_*(lambda_Adot_HR_iARD[1]*lambda_Adot_HR_iARD[1]
                                        + 2.0*lambda_Adot_HR_iARD[2]*lambda_Adot_HR_iARD[0])
                            ); 
                        
        
    LambdaDot_IOK.zz() = alpha_*(
                                lambda_Adot_HR_iARD[2]
                                -beta_*(lambda_Adot_HR_iARD[2]*lambda_Adot_HR_iARD[2]
                                        + 2.0*lambda_Adot_HR_iARD[0]*lambda_Adot_HR_iARD[1])
                            ); 
}
