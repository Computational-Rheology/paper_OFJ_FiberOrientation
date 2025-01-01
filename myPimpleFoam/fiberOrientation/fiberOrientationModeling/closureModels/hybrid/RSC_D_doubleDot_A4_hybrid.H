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

Literature
   S. G. Advani and C. L. Tucker,
   “Closure approximations for three‐dimensional structure tensors,” 
   J. Rheol. (N. Y. N. Y)., vol. 34, no. 3, pp. 367–386, 1990, doi: 10.1122/1.550133.

   J. H. Phelps and C. L. Tucker, 
   “An anisotropic rotary diffusion model for fiber orientation in short- and long-fiber thermoplastics,”
    J. Nonnewton. Fluid Mech., vol. 156, no. 3, pp. 165–176, 2009, doi: 10.1016/j.jnnfm.2008.08.002.

\*----------------------------------------------------------------------------*/


void Foam::fiberOrientation::closureModels::hybrid::RSC_D_doubleDot_A4_hybrid
(
   Foam::symmTensor& result,
   const Foam::symmTensor& D,
   const Foam::symmTensor& A,
   const Foam::vector& eigenValue,
   const Foam::tensor& eigenVector,
   const Foam::scalar& k
)
{
    const scalar tmp0 = 1.0 - k;
    const scalar tmp1 = eigenVector.xx()*eigenVector.xx()*eigenVector.xx()*eigenVector.xx();
    const scalar tmp2 = eigenVector.yx()*eigenVector.yx()*eigenVector.yx()*eigenVector.yx();
    const scalar tmp3 = eigenVector.zx()*eigenVector.zx()*eigenVector.zx()*eigenVector.zx();
    const scalar tmp4 = tmp1 + tmp2 + tmp3;
    const scalar tmp5 = A.yz()*A.yz();
    const scalar tmp6 = A.xy()*A.xy();
    const scalar tmp7 = 27.0*A.zz();
    const scalar tmp8 = A.xz()*A.xz();
    const scalar tmp9 = A.xz()*A.yz();
    const scalar tmp10 = 27.0*tmp5*A.xx() + tmp6*tmp7 - tmp7*A.xx()*A.yy() + 27.0*tmp8*A.yy() - 54.0*tmp9*A.xy();
    const scalar tmp11 = tmp10 + 1.0;
    const scalar tmp12 = -tmp10;
    const scalar tmp13 = tmp11*(A.xx()*A.xx()) + tmp12*(0.85714285714285714*A.xx() - 0.085714285714285714);
    const scalar tmp14 = eigenVector.xx()*eigenVector.xx()*eigenVector.xx();
    const scalar tmp15 = eigenVector.yx()*eigenVector.yx()*eigenVector.yx();
    const scalar tmp16 = eigenVector.zx()*eigenVector.zx()*eigenVector.zx();
    const scalar tmp17 = tmp14*eigenVector.xy() + tmp15*eigenVector.yy() + tmp16*eigenVector.zy();
    const scalar tmp18 = 0.42857142857142857*tmp12;
    const scalar tmp19 = tmp18*A.xy();
    const scalar tmp20 = tmp11*A.xx();
    const scalar tmp21 = tmp19 + tmp20*A.xy();
    const scalar tmp22 = tmp17*tmp21;
    const scalar tmp23 = tmp14*eigenVector.xz() + tmp15*eigenVector.yz() + tmp16*eigenVector.zz();
    const scalar tmp24 = tmp18*A.xz();
    const scalar tmp25 = tmp20*A.xz() + tmp24;
    const scalar tmp26 = tmp23*tmp25;
    const scalar tmp27 = eigenVector.xx()*eigenVector.xx();
    const scalar tmp28 = eigenVector.xy()*eigenVector.xy();
    const scalar tmp29 = eigenVector.yx()*eigenVector.yx();
    const scalar tmp30 = eigenVector.yy()*eigenVector.yy();
    const scalar tmp31 = eigenVector.zx()*eigenVector.zx();
    const scalar tmp32 = eigenVector.zy()*eigenVector.zy();
    const scalar tmp33 = tmp27*tmp28 + tmp29*tmp30 + tmp31*tmp32;
    const scalar tmp34 = 0.14285714285714286*A.yy();
    const scalar tmp35 = 0.14285714285714286*A.xx() - 0.028571428571428571;
    const scalar tmp36 = tmp12*(tmp34 + tmp35);
    const scalar tmp37 = tmp20*A.yy() + tmp36;
    const scalar tmp38 = tmp33*tmp37;
    const scalar tmp39 = eigenVector.xz()*eigenVector.xz();
    const scalar tmp40 = eigenVector.yz()*eigenVector.yz();
    const scalar tmp41 = eigenVector.zz()*eigenVector.zz();
    const scalar tmp42 = tmp27*tmp39 + tmp29*tmp40 + tmp31*tmp41;
    const scalar tmp43 = 0.14285714285714286*A.zz();
    const scalar tmp44 = tmp12*(tmp35 + tmp43);
    const scalar tmp45 = tmp20*A.zz() + tmp44;
    const scalar tmp46 = tmp42*tmp45;
    const scalar tmp47 = tmp27*eigenVector.xy()*eigenVector.xz() + tmp29*eigenVector.yy()*eigenVector.yz() + tmp31*eigenVector.zy()*eigenVector.zz();
    const scalar tmp48 = 0.14285714285714286*tmp12;
    const scalar tmp49 = tmp48*A.yz();
    const scalar tmp50 = tmp20*A.yz() + tmp49;
    const scalar tmp51 = tmp47*tmp50;
    const scalar tmp52 = tmp48*A.xy();
    const scalar tmp53 = tmp11*A.zz();
    const scalar tmp54 = tmp52 + tmp53*A.xy();
    const scalar tmp55 = tmp11*A.yy();
    const scalar tmp56 = tmp19 + tmp55*A.xy();
    const scalar tmp57 = tmp33*tmp56;
    const scalar tmp58 = tmp11*tmp6 + tmp36;
    const scalar tmp59 = 2.0*tmp17;
    const scalar tmp60 = tmp48*A.xz();
    const scalar tmp61 = tmp11*A.xy();
    const scalar tmp62 = tmp60 + tmp61*A.yz();
    const scalar tmp63 = 2.0*tmp47;
    const scalar tmp64 = tmp62*tmp63;
    const scalar tmp65 = tmp49 + tmp61*A.xz();
    const scalar tmp66 = 2.0*tmp23;
    const scalar tmp67 = -tmp14*eigenValue.x()*eigenVector.xy() - tmp15*eigenValue.y()*eigenVector.yy() - tmp16*eigenValue.z()*eigenVector.zy();
    const scalar tmp68 = 2.0*D.xy();
    const scalar tmp69 = tmp55*A.xz() + tmp60;
    const scalar tmp70 = tmp24 + tmp53*A.xz();
    const scalar tmp71 = tmp42*tmp70;
    const scalar tmp72 = tmp11*tmp8 + tmp44;
    const scalar tmp73 = tmp11*tmp9 + tmp52;
    const scalar tmp74 = tmp63*tmp73;
    const scalar tmp75 = -tmp14*eigenValue.x()*eigenVector.xz() - tmp15*eigenValue.y()*eigenVector.yz() - tmp16*eigenValue.z()*eigenVector.zz();
    const scalar tmp76 = 2.0*D.xz();
    const scalar tmp77 = tmp18*A.yz();
    const scalar tmp78 = tmp55*A.yz() + tmp77;
    const scalar tmp79 = tmp53*A.yz() + tmp77;
    const scalar tmp80 = tmp12*(tmp34 + tmp43 - 0.028571428571428571);
    const scalar tmp81 = tmp11*tmp5 + tmp80;
    const scalar tmp82 = -tmp27*eigenValue.x()*eigenVector.xy()*eigenVector.xz() - tmp29*eigenValue.y()*eigenVector.yy()*eigenVector.yz() - tmp31*eigenValue.z()*eigenVector.zy()*eigenVector.zz();
    const scalar tmp83 = 2.0*D.yz();
    const scalar tmp84 = tmp11*(A.yy()*A.yy()) + tmp12*(0.85714285714285714*A.yy() - 0.085714285714285714);
    const scalar tmp85 = tmp55*A.zz() + tmp80;
    const scalar tmp86 = -tmp27*tmp28*eigenValue.x() - tmp29*tmp30*eigenValue.y() - tmp31*tmp32*eigenValue.z();
    const scalar tmp87 = tmp11*(A.zz()*A.zz()) + tmp12*(0.85714285714285714*A.zz() - 0.085714285714285714);
    const scalar tmp88 = -tmp27*tmp39*eigenValue.x() - tmp29*tmp40*eigenValue.y() - tmp31*tmp41*eigenValue.z();
    const scalar tmp89 = tmp39*eigenVector.xx()*eigenVector.xy() + tmp40*eigenVector.yx()*eigenVector.yy() + tmp41*eigenVector.zx()*eigenVector.zy();
    const scalar tmp90 = tmp54*tmp89;
    const scalar tmp91 = eigenVector.xy()*eigenVector.xy()*eigenVector.xy();
    const scalar tmp92 = eigenVector.yy()*eigenVector.yy()*eigenVector.yy();
    const scalar tmp93 = eigenVector.zy()*eigenVector.zy()*eigenVector.zy();
    const scalar tmp94 = tmp91*eigenVector.xx() + tmp92*eigenVector.yx() + tmp93*eigenVector.zx();
    const scalar tmp95 = tmp56*tmp94;
    const scalar tmp96 = 2.0*tmp33;
    const scalar tmp97 = tmp28*eigenVector.xx()*eigenVector.xz() + tmp30*eigenVector.yx()*eigenVector.yz() + tmp32*eigenVector.zx()*eigenVector.zz();
    const scalar tmp98 = 2.0*tmp97;
    const scalar tmp99 = tmp62*tmp98;
    const scalar tmp100 = tmp63*tmp65;
    const scalar tmp101 = tmp21*tmp33;
    const scalar tmp102 = tmp25*tmp47;
    const scalar tmp103 = tmp78*tmp97;
    const scalar tmp104 = -tmp91*eigenValue.x()*eigenVector.xx() - tmp92*eigenValue.y()*eigenVector.yx() - tmp93*eigenValue.z()*eigenVector.zx();
    const scalar tmp105 = tmp70*tmp89;
    const scalar tmp106 = tmp73*tmp98;
    const scalar tmp107 = tmp79*tmp89;
    const scalar tmp108 = -tmp28*eigenValue.x()*eigenVector.xx()*eigenVector.xz() - tmp30*eigenValue.y()*eigenVector.yx()*eigenVector.yz() - tmp32*eigenValue.z()*eigenVector.zx()*eigenVector.zz();
    const scalar tmp109 = -tmp39*eigenValue.x()*eigenVector.xx()*eigenVector.xy() - tmp40*eigenValue.y()*eigenVector.yx()*eigenVector.yy() - tmp41*eigenValue.z()*eigenVector.zx()*eigenVector.zy();
    const scalar tmp110 = tmp68*(tmp0*(-tmp100 - tmp22 - tmp58*tmp96 - tmp86 - tmp90 - tmp95 - tmp99) + tmp58) + tmp76*(tmp0*(-tmp105 - tmp106 - tmp17*tmp25 - tmp63*tmp72 - tmp65*tmp96 - tmp69*tmp94 - tmp82) + tmp65) + tmp83*(tmp0*(-tmp107 - tmp108 - tmp17*tmp50 - tmp62*tmp96 - tmp74 - tmp78*tmp94 - tmp81*tmp98) + tmp62) + (tmp0*(-2.0*tmp101 - 2.0*tmp102 - tmp13*tmp17 - tmp37*tmp94 - tmp45*tmp89 - tmp50*tmp98 - tmp67) + tmp21)*D.xx() + (tmp0*(-2.0*tmp103 - tmp104 - tmp17*tmp37 - 2.0*tmp57 - tmp63*tmp69 - tmp84*tmp94 - tmp85*tmp89) + tmp56)*D.yy() + (tmp0*(-tmp109 - tmp17*tmp45 - tmp54*tmp96 - tmp63*tmp70 - tmp79*tmp98 - tmp85*tmp94 - tmp87*tmp89) + tmp54)*D.zz();
    const scalar tmp111 = tmp69*tmp97;
    const scalar tmp112 = eigenVector.xz()*eigenVector.xz()*eigenVector.xz();
    const scalar tmp113 = eigenVector.yz()*eigenVector.yz()*eigenVector.yz();
    const scalar tmp114 = eigenVector.zz()*eigenVector.zz()*eigenVector.zz();
    const scalar tmp115 = tmp112*eigenVector.xx() + tmp113*eigenVector.yx() + tmp114*eigenVector.zx();
    const scalar tmp116 = tmp115*tmp70;
    const scalar tmp117 = 2.0*tmp42;
    const scalar tmp118 = 2.0*tmp89;
    const scalar tmp119 = tmp118*tmp73;
    const scalar tmp120 = tmp21*tmp47;
    const scalar tmp121 = tmp25*tmp42;
    const scalar tmp122 = -tmp112*eigenValue.x()*eigenVector.xx() - tmp113*eigenValue.y()*eigenVector.yx() - tmp114*eigenValue.z()*eigenVector.zx();
    const scalar tmp123 = tmp56*tmp97;
    const scalar tmp124 = tmp118*tmp62;
    const scalar tmp125 = tmp68*(tmp0*(-tmp115*tmp54 - tmp117*tmp65 - tmp123 - tmp124 - tmp21*tmp23 - tmp58*tmp63 - tmp82) + tmp65) + tmp76*(tmp0*(-tmp100 - tmp111 - tmp116 - tmp117*tmp72 - tmp119 - tmp26 - tmp88) + tmp72) + tmp83*(tmp0*(-tmp103 - tmp109 - tmp115*tmp79 - tmp117*tmp73 - tmp118*tmp81 - tmp23*tmp50 - tmp64) + tmp73) + (tmp0*(-2.0*tmp107 - tmp115*tmp87 - tmp122 - tmp23*tmp45 - tmp54*tmp63 - 2.0*tmp71 - tmp85*tmp97) + tmp70)*D.zz() + (tmp0*(-tmp108 - tmp115*tmp85 - tmp117*tmp69 - tmp118*tmp78 - tmp23*tmp37 - tmp56*tmp63 - tmp84*tmp97) + tmp69)*D.yy() + (tmp0*(-tmp115*tmp45 - tmp118*tmp50 - 2.0*tmp120 - 2.0*tmp121 - tmp13*tmp23 - tmp37*tmp97 - tmp75) + tmp25)*D.xx();
    const scalar tmp126 = eigenVector.xy()*eigenVector.xy()*eigenVector.xy()*eigenVector.xy();
    const scalar tmp127 = eigenVector.yy()*eigenVector.yy()*eigenVector.yy()*eigenVector.yy();
    const scalar tmp128 = eigenVector.zy()*eigenVector.zy()*eigenVector.zy()*eigenVector.zy();
    const scalar tmp129 = tmp126 + tmp127 + tmp128;
    const scalar tmp130 = tmp91*eigenVector.xz() + tmp92*eigenVector.yz() + tmp93*eigenVector.zz();
    const scalar tmp131 = tmp130*tmp78;
    const scalar tmp132 = tmp28*tmp39 + tmp30*tmp40 + tmp32*tmp41;
    const scalar tmp133 = tmp132*tmp85;
    const scalar tmp134 = 2.0*tmp94;
    const scalar tmp135 = 2.0*tmp130;
    const scalar tmp136 = tmp65*tmp98;
    const scalar tmp137 = tmp132*tmp79;
    const scalar tmp138 = -tmp91*eigenValue.x()*eigenVector.xz() - tmp92*eigenValue.y()*eigenVector.yz() - tmp93*eigenValue.z()*eigenVector.zz();
    const scalar tmp139 = -tmp28*tmp39*eigenValue.x() - tmp30*tmp40*eigenValue.y() - tmp32*tmp41*eigenValue.z();
    const scalar tmp140 = tmp112*eigenVector.xy() + tmp113*eigenVector.yy() + tmp114*eigenVector.zy();
    const scalar tmp141 = tmp140*tmp79;
    const scalar tmp142 = 2.0*tmp132;
    const scalar tmp143 = tmp132*tmp78;
    const scalar tmp144 = -tmp112*eigenValue.x()*eigenVector.xy() - tmp113*eigenValue.y()*eigenVector.yy() - tmp114*eigenValue.z()*eigenVector.zy();
    const scalar tmp145 = tmp118*tmp65;
    const scalar tmp146 = tmp68*(tmp0*(-tmp108 - tmp120 - tmp130*tmp56 - tmp140*tmp54 - tmp142*tmp62 - tmp145 - tmp58*tmp98) + tmp62) + tmp76*(tmp0*(-tmp102 - tmp109 - tmp118*tmp72 - tmp130*tmp69 - tmp136 - tmp140*tmp70 - tmp142*tmp73) + tmp73) + tmp83*(tmp0*(-tmp119 - tmp131 - tmp139 - tmp141 - tmp142*tmp81 - tmp51 - tmp99) + tmp81) + (tmp0*(-2.0*tmp105 - tmp130*tmp85 - 2.0*tmp137 - tmp140*tmp87 - tmp144 - tmp45*tmp47 - tmp54*tmp98) + tmp79)*D.zz() + (tmp0*(-tmp118*tmp25 - tmp13*tmp47 - tmp130*tmp37 - tmp140*tmp45 - tmp142*tmp50 - tmp21*tmp98 - tmp82) + tmp50)*D.xx() + (tmp0*(-tmp118*tmp69 - 2.0*tmp123 - tmp130*tmp84 - tmp138 - tmp140*tmp85 - 2.0*tmp143 - tmp37*tmp47) + tmp78)*D.yy();
    const scalar tmp147 = eigenVector.xz()*eigenVector.xz()*eigenVector.xz()*eigenVector.xz();
    const scalar tmp148 = eigenVector.yz()*eigenVector.yz()*eigenVector.yz()*eigenVector.yz();
    const scalar tmp149 = eigenVector.zz()*eigenVector.zz()*eigenVector.zz()*eigenVector.zz();
    const scalar tmp150 = tmp147 + tmp148 + tmp149;
    const scalar tmp151 = 2.0*tmp115;
    const scalar tmp152 = 2.0*tmp140;

    result.xx() = tmp68*(tmp0*(-tmp21*tmp4 - tmp42*tmp54 - tmp57 - tmp58*tmp59 - tmp64 - tmp65*tmp66 - tmp67) + tmp21) + tmp76*(tmp0*(-tmp25*tmp4 - tmp33*tmp69 - tmp59*tmp65 - tmp66*tmp72 - tmp71 - tmp74 - tmp75) + tmp25) + tmp83*(tmp0*(-tmp33*tmp78 - tmp4*tmp50 - tmp42*tmp79 - tmp59*tmp62 - tmp63*tmp81 - tmp66*tmp73 - tmp82) + tmp50) + (tmp0*(-tmp33*tmp84 - tmp37*tmp4 - tmp42*tmp85 - tmp56*tmp59 - tmp63*tmp78 - tmp66*tmp69 - tmp86) + tmp37)*D.yy() + (tmp0*(-tmp33*tmp85 - tmp4*tmp45 - tmp42*tmp87 - tmp54*tmp59 - tmp63*tmp79 - tmp66*tmp70 - tmp88) + tmp45)*D.zz() + (tmp0*(tmp1*eigenValue.x() - tmp13*tmp4 + tmp2*eigenValue.y() - 2.0*tmp22 - 2.0*tmp26 + tmp3*eigenValue.z() - tmp38 - tmp46 - 2.0*tmp51) + tmp13)*D.xx();
    result.xy() = tmp110;
    result.xz() = tmp125;
    result.yy() = tmp68*(tmp0*(-tmp101 - tmp104 - tmp129*tmp56 - tmp132*tmp54 - tmp134*tmp58 - tmp135*tmp62 - tmp136) + tmp56) + tmp76*(tmp0*(-tmp108 - tmp129*tmp69 - tmp132*tmp70 - tmp134*tmp65 - tmp135*tmp73 - tmp25*tmp33 - tmp72*tmp98) + tmp69) + tmp83*(tmp0*(-tmp106 - tmp129*tmp78 - tmp134*tmp62 - tmp135*tmp81 - tmp137 - tmp138 - tmp33*tmp50) + tmp78) + (tmp0*(-tmp129*tmp37 - tmp13*tmp33 - tmp132*tmp45 - tmp134*tmp21 - tmp135*tmp50 - tmp25*tmp98 - tmp86) + tmp37)*D.xx() + (tmp0*(-tmp129*tmp85 - tmp132*tmp87 - tmp134*tmp54 - tmp135*tmp79 - tmp139 - tmp33*tmp45 - tmp70*tmp98) + tmp85)*D.zz() + (tmp0*(-2.0*tmp111 + tmp126*eigenValue.x() + tmp127*eigenValue.y() + tmp128*eigenValue.z() - tmp129*tmp84 - 2.0*tmp131 - tmp133 - tmp38 - 2.0*tmp95) + tmp84)*D.yy();
    result.yz() = tmp146;
    result.zz() = tmp68*(tmp0*(-tmp109 - tmp118*tmp58 - tmp132*tmp56 - tmp150*tmp54 - tmp151*tmp65 - tmp152*tmp62 - tmp21*tmp42) + tmp54) + tmp76*(tmp0*(-tmp121 - tmp122 - tmp132*tmp69 - tmp145 - tmp150*tmp70 - tmp151*tmp72 - tmp152*tmp73) + tmp70) + tmp83*(tmp0*(-tmp124 - tmp143 - tmp144 - tmp150*tmp79 - tmp151*tmp73 - tmp152*tmp81 - tmp42*tmp50) + tmp79) + (tmp0*(-tmp118*tmp21 - tmp13*tmp42 - tmp132*tmp37 - tmp150*tmp45 - tmp151*tmp25 - tmp152*tmp50 - tmp88) + tmp45)*D.xx() + (tmp0*(-tmp118*tmp56 - tmp132*tmp84 - tmp139 - tmp150*tmp85 - tmp151*tmp69 - tmp152*tmp78 - tmp37*tmp42) + tmp85)*D.yy() + (tmp0*(-2.0*tmp116 - tmp133 - 2.0*tmp141 + tmp147*eigenValue.x() + tmp148*eigenValue.y() + tmp149*eigenValue.z() - tmp150*tmp87 - tmp46 - 2.0*tmp90) + tmp87)*D.zz();
}
