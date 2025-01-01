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


Description
    The following function was generated with symbolic computation
\*----------------------------------------------------------------------------*/

void Foam::fiberOrientation::closureModels::hybrid::D_doubleDot_A4_hybrid
(
    Foam::symmTensor& result, 
    const Foam::symmTensor& D, 
    const Foam::symmTensor& A
)
{
    const scalar tmp0 = A.yz()*A.yz();
    const scalar tmp1 = A.xy()*A.xy();
    const scalar tmp2 = 27.0*A.zz();
    const scalar tmp3 = A.xz()*A.xz();
    const scalar tmp4 = 27.0*tmp0*A.xx() + tmp1*tmp2 - tmp2*A.xx()*A.yy() + 27.0*tmp3*A.yy() - 54.0*A.xy()*A.xz()*A.yz();
    const scalar tmp5 = tmp4 + 1.0;
    const scalar tmp6 = -tmp4;
    const scalar tmp7 = 0.42857142857142857*tmp6;
    const scalar tmp8 = tmp7*A.xy();
    const scalar tmp9 = tmp5*A.xx();
    const scalar tmp10 = tmp8 + tmp9*A.xy();
    const scalar tmp11 = 2.0*D.xy();
    const scalar tmp12 = tmp7*A.xz();
    const scalar tmp13 = tmp12 + tmp9*A.xz();
    const scalar tmp14 = 2.0*D.xz();
    const scalar tmp15 = 0.14285714285714286*tmp6;
    const scalar tmp16 = tmp15*A.yz();
    const scalar tmp17 = tmp16 + tmp9*A.yz();
    const scalar tmp18 = 2.0*D.yz();
    const scalar tmp19 = 0.14285714285714286*A.yy();
    const scalar tmp20 = 0.14285714285714286*A.xx() - 0.028571428571428571;
    const scalar tmp21 = tmp6*(tmp19 + tmp20);
    const scalar tmp22 = tmp21 + tmp9*A.yy();
    const scalar tmp23 = 0.14285714285714286*A.zz();
    const scalar tmp24 = tmp6*(tmp20 + tmp23);
    const scalar tmp25 = tmp24 + tmp9*A.zz();
    const scalar tmp26 = tmp15*A.xy();
    const scalar tmp27 = tmp5*A.xy();
    const scalar tmp28 = tmp26 + tmp27*A.zz();
    const scalar tmp29 = tmp27*A.yy() + tmp8;
    const scalar tmp30 = tmp15*A.xz();
    const scalar tmp31 = tmp27*A.yz() + tmp30;
    const scalar tmp32 = tmp16 + tmp27*A.xz();
    const scalar tmp33 = tmp10*D.xx() + tmp11*(tmp1*tmp5 + tmp21) + tmp14*tmp32 + tmp18*tmp31 + tmp28*D.zz() + tmp29*D.yy();
    const scalar tmp34 = tmp5*A.xz();
    const scalar tmp35 = tmp30 + tmp34*A.yy();
    const scalar tmp36 = tmp12 + tmp34*A.zz();
    const scalar tmp37 = tmp26 + tmp34*A.yz();
    const scalar tmp38 = tmp11*tmp32 + tmp13*D.xx() + tmp14*(tmp24 + tmp3*tmp5) + tmp18*tmp37 + tmp35*D.yy() + tmp36*D.zz();
    const scalar tmp39 = tmp7*A.yz();
    const scalar tmp40 = tmp5*A.yy();
    const scalar tmp41 = tmp39 + tmp40*A.yz();
    const scalar tmp42 = tmp6*(tmp19 + tmp23 - 0.028571428571428571);
    const scalar tmp43 = tmp40*A.zz() + tmp42;
    const scalar tmp44 = tmp39 + tmp5*A.yz()*A.zz();
    const scalar tmp45 = tmp11*tmp31 + tmp14*tmp37 + tmp17*D.xx() + tmp18*(tmp0*tmp5 + tmp42) + tmp41*D.yy() + tmp44*D.zz();

    result.xx() = tmp10*tmp11 + tmp13*tmp14 + tmp17*tmp18 + tmp22*D.yy() + tmp25*D.zz() + D.xx()*(tmp5*(A.xx()*A.xx()) + tmp6*(0.85714285714285714*A.xx() - 0.085714285714285714));
    result.xy() = tmp33;
    result.xz() = tmp38;
    result.yy() = tmp11*tmp29 + tmp14*tmp35 + tmp18*tmp41 + tmp22*D.xx() + tmp43*D.zz() + D.yy()*(tmp5*(A.yy()*A.yy()) + tmp6*(0.85714285714285714*A.yy() - 0.085714285714285714));
    result.yz() = tmp45;
    result.zz() = tmp11*tmp28 + tmp14*tmp36 + tmp18*tmp44 + tmp25*D.xx() + tmp43*D.yy() + D.zz()*(tmp5*(A.zz()*A.zz()) + tmp6*(0.85714285714285714*A.zz() - 0.085714285714285714));
}
