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
    B. E. VerWeyst.
    "Numerical predictions of flow-induced fiber orientation in 3-D geometries". PhD thesis,
    University of Illinois at Urbana-Champaign, Urbana, IL, 1998.
    
\*-------------------------------------------------------------------------------------------------------*/

void Foam::fiberOrientation::closureModels::ORE::RSC_D_doubleDot_A4_ORE
(
    Foam::symmTensor& result,
    const Foam::symmTensor& D,
    const Foam::vector& eigenValue,
    const Foam::tensor& eigenVector,
    const Foam::scalar&  k
)
{
    const scalar tmp0 = 1.0 - k;
    const scalar tmp1 = eigenVector.xx()*eigenVector.xx()*eigenVector.xx()*eigenVector.xx();
    const scalar tmp2 = eigenVector.yx()*eigenVector.yx()*eigenVector.yx()*eigenVector.yx();
    const scalar tmp3 = eigenVector.zx()*eigenVector.zx()*eigenVector.zx()*eigenVector.zx();
    const scalar tmp4 = tmp1 + tmp2 + tmp3;
    const scalar tmp5 = eigenValue.y()*eigenValue.y();
    const scalar tmp6 = eigenValue.y()*eigenValue.y()*eigenValue.y();
    const scalar tmp7 = eigenValue.y()*eigenValue.y()*eigenValue.y()*eigenValue.y();
    const scalar tmp8 = eigenValue.z()*eigenValue.z();
    const scalar tmp9 = eigenValue.z()*eigenValue.z()*eigenValue.z();
    const scalar tmp10 = eigenValue.z()*eigenValue.z()*eigenValue.z()*eigenValue.z();
    const scalar tmp11 = eigenValue.y()*eigenValue.z();
    const scalar tmp12 = tmp8*eigenValue.y();
    const scalar tmp13 = tmp9*eigenValue.y();
    const scalar tmp14 = tmp5*eigenValue.z();
    const scalar tmp15 = tmp6*eigenValue.z();
    const scalar tmp16 = tmp5*tmp8;
    const scalar tmp17 = 2.28476531637958*tmp10 + 11.8273285968852*tmp11 - 15.9120667157641*tmp12 + 7.74683751713295*tmp13 - 15.1515872606307*tmp14 + 7.4814687062444101*tmp15 + 9.3252034345266104*tmp16 + 8.4367774677832497*tmp5 - 8.6389141928401596*tmp6 + 3.59772251134254*tmp7 + 6.8815395205804402*tmp8 - 6.4872893364192601*tmp9 - 3.0370993982540599*eigenValue.y() - 3.3152722974214601*eigenValue.z() + 0.63625679688068704;
    const scalar tmp18 = 0.50804138736663695*tmp10 + 11.958956233232*tmp11 - 10.9582626069691*tmp12 + 3.4790151056743901*tmp13 - 20.727799468413199*tmp14 + 11.749291117702599*tmp15 + 9.8159838971674809*tmp16 + 11.342092427815899*tmp5 - 12.3875632855619*tmp6 + 4.8836659777148901*tmp7 + 3.84459692420086*tmp8 - 2.1162321447100401*tmp9 - 4.4797087319373796*eigenValue.y() - 1.87266296373814*eigenValue.z() + 0.63625679688068704;
    const scalar tmp19 = 3.43213840334779*tmp10 + 34.319901891698699*tmp11 - 37.702911802938402*tmp12 + 15.2650686148651*tmp13 - 50.275643192748497*tmp14 + 26.113491400537502*tmp15 + 27.3346798054488*tmp16 + 25.868475525388401*tmp5 - 26.963691523971601*tmp6 + 10.611741806606*tmp7 + 13.829469912194*tmp8 - 10.8801761133174*tmp9 - 12.2570587036254*eigenValue.y() - 9.1219650978269193*eigenValue.z() + 2.7405328956025299;
    const scalar tmp20 = eigenVector.xx()*eigenVector.xx();
    const scalar tmp21 = -0.82770723716742345*tmp10 - 17.225764764022749*tmp11 - 5.4986231017032701*tmp13 - 15.190656905997844*tmp15 - 13.912730134044835*tmp16 + 27.925927700265497*tmp5*eigenValue.z() - 14.386895242710526*tmp5 + 15.35617030834667*tmp6 - 5.9488426364891751*tmp7 + 16.374553847071702*tmp8*eigenValue.y() - 5.3962636579072099*tmp8 + 3.2545594608040904*tmp9 + 5.8498340186543603*eigenValue.y() + 3.8396778820717996*eigenValue.z() - 0.87026644780126494;
    const scalar tmp22 = tmp20*tmp21;
    const scalar tmp23 = eigenVector.zx()*eigenVector.zx();
    const scalar tmp24 = 6.0*tmp23;
    const scalar tmp25 = eigenVector.yx()*eigenVector.yx();
    const scalar tmp26 = -2.6044311661803663*tmp10 - 17.09413712767595*tmp11 - 9.7664455131618304*tmp13 - 10.922834494539657*tmp15 - 13.421949671403965*tmp16 + 22.349715492483*tmp5*eigenValue.z() - 11.481580282677875*tmp5 + 11.607521215624931*tmp6 - 4.662899170116825*tmp7 + 21.3283579558667*tmp8*eigenValue.y() - 8.4332062542867909*tmp8 + 7.62561665251331*tmp9 + 5.4072246849710401*eigenValue.y() + 4.2822872157551197*eigenValue.z() - 0.87026644780126494;
    const scalar tmp27 = tmp20*tmp26;
    const scalar tmp28 = 0.31966584980078649*tmp10 + 5.2668085307907493*tmp11 - 5.416291240102602*tmp12 + 2.0196079960288804*tmp13 - 7.1981282318522979*tmp14 + 3.4413657882952453*tmp15 + 4.0967462368773546*tmp16 + 3.0448028148946253*tmp5 - 2.9686070227847718*tmp6 + 1.065176658774285*tmp7 + 1.5516667337063499*tmp8 - 1.1383273160940499*tmp9 - 1.3701252867169806*eigenValue.y() - 0.9670149183336596*eigenValue.z() + 0.23400965092057779;
    const scalar tmp29 = tmp25*tmp28;
    const scalar tmp30 = tmp1*tmp19 + tmp17*tmp2 + tmp18*tmp3 + tmp22*tmp24 + tmp24*tmp29 + 6.0*tmp25*tmp27;
    const scalar tmp31 = eigenVector.xx()*eigenVector.xx()*eigenVector.xx();
    const scalar tmp32 = tmp31*eigenVector.xy();
    const scalar tmp33 = eigenVector.yx()*eigenVector.yx()*eigenVector.yx();
    const scalar tmp34 = tmp33*eigenVector.yy();
    const scalar tmp35 = eigenVector.zx()*eigenVector.zx()*eigenVector.zx();
    const scalar tmp36 = tmp35*eigenVector.zy();
    const scalar tmp37 = tmp32 + tmp34 + tmp36;
    const scalar tmp38 = tmp21*tmp23;
    const scalar tmp39 = eigenVector.xx()*eigenVector.xy();
    const scalar tmp40 = 3.0*tmp39;
    const scalar tmp41 = eigenVector.zx()*eigenVector.zy();
    const scalar tmp42 = 3.0*tmp41;
    const scalar tmp43 = tmp25*tmp26;
    const scalar tmp44 = eigenVector.yx()*eigenVector.yy();
    const scalar tmp45 = 3.0*tmp44;
    const scalar tmp46 = tmp23*tmp28;
    const scalar tmp47 = tmp17*tmp34 + tmp18*tmp36 + tmp19*tmp32 + tmp22*tmp42 + tmp27*tmp45 + tmp29*tmp42 + tmp38*tmp40 + tmp40*tmp43 + tmp45*tmp46;
    const scalar tmp48 = tmp37*tmp47;
    const scalar tmp49 = tmp31*eigenVector.xz();
    const scalar tmp50 = tmp33*eigenVector.yz();
    const scalar tmp51 = tmp35*eigenVector.zz();
    const scalar tmp52 = tmp49 + tmp50 + tmp51;
    const scalar tmp53 = eigenVector.xx()*eigenVector.xz();
    const scalar tmp54 = 3.0*tmp53;
    const scalar tmp55 = eigenVector.zx()*eigenVector.zz();
    const scalar tmp56 = 3.0*tmp55;
    const scalar tmp57 = eigenVector.yx()*eigenVector.yz();
    const scalar tmp58 = 3.0*tmp57;
    const scalar tmp59 = tmp17*tmp50 + tmp18*tmp51 + tmp19*tmp49 + tmp22*tmp56 + tmp27*tmp58 + tmp29*tmp56 + tmp38*tmp54 + tmp43*tmp54 + tmp46*tmp58;
    const scalar tmp60 = tmp52*tmp59;
    const scalar tmp61 = eigenVector.xy()*eigenVector.xy();
    const scalar tmp62 = tmp20*tmp61;
    const scalar tmp63 = eigenVector.yy()*eigenVector.yy();
    const scalar tmp64 = tmp25*tmp63;
    const scalar tmp65 = eigenVector.zy()*eigenVector.zy();
    const scalar tmp66 = tmp23*tmp65;
    const scalar tmp67 = tmp62 + tmp64 + tmp66;
    const scalar tmp68 = 4.0*tmp41;
    const scalar tmp69 = tmp21*tmp39;
    const scalar tmp70 = 4.0*tmp26;
    const scalar tmp71 = tmp28*tmp44;
    const scalar tmp72 = tmp17*tmp64 + tmp18*tmp66 + tmp19*tmp62 + tmp22*tmp65 + tmp27*tmp63 + tmp29*tmp65 + tmp38*tmp61 + tmp39*tmp44*tmp70 + tmp43*tmp61 + tmp46*tmp63 + tmp68*tmp69 + tmp68*tmp71;
    const scalar tmp73 = tmp67*tmp72;
    const scalar tmp74 = eigenVector.xz()*eigenVector.xz();
    const scalar tmp75 = tmp20*tmp74;
    const scalar tmp76 = eigenVector.yz()*eigenVector.yz();
    const scalar tmp77 = tmp25*tmp76;
    const scalar tmp78 = eigenVector.zz()*eigenVector.zz();
    const scalar tmp79 = tmp23*tmp78;
    const scalar tmp80 = tmp75 + tmp77 + tmp79;
    const scalar tmp81 = 4.0*tmp55;
    const scalar tmp82 = tmp21*tmp53;
    const scalar tmp83 = tmp28*tmp57;
    const scalar tmp84 = tmp17*tmp77 + tmp18*tmp79 + tmp19*tmp75 + tmp22*tmp78 + tmp27*tmp76 + tmp29*tmp78 + tmp38*tmp74 + tmp43*tmp74 + tmp46*tmp76 + tmp53*tmp57*tmp70 + tmp81*tmp82 + tmp81*tmp83;
    const scalar tmp85 = tmp80*tmp84;
    const scalar tmp86 = eigenVector.xy()*eigenVector.xz();
    const scalar tmp87 = tmp20*tmp86;
    const scalar tmp88 = eigenVector.yy()*eigenVector.yz();
    const scalar tmp89 = tmp25*tmp88;
    const scalar tmp90 = eigenVector.zy()*eigenVector.zz();
    const scalar tmp91 = tmp23*tmp90;
    const scalar tmp92 = tmp87 + tmp89 + tmp91;
    const scalar tmp93 = 2.0*tmp55;
    const scalar tmp94 = 2.0*tmp41;
    const scalar tmp95 = 2.0*tmp26;
    const scalar tmp96 = tmp39*tmp95;
    const scalar tmp97 = tmp44*tmp95;
    const scalar tmp98 = tmp17*tmp89 + tmp18*tmp91 + tmp19*tmp87 + tmp22*tmp90 + tmp27*tmp88 + tmp29*tmp90 + tmp38*tmp86 + tmp43*tmp86 + tmp46*tmp88 + tmp53*tmp97 + tmp57*tmp96 + tmp69*tmp93 + tmp71*tmp93 + tmp82*tmp94 + tmp83*tmp94;
    const scalar tmp99 = tmp92*tmp98;
    const scalar tmp100 = 2.0*tmp99;
    const scalar tmp101 = eigenVector.yy()*eigenVector.yy()*eigenVector.yy()*eigenVector.yy();
    const scalar tmp102 = eigenVector.zy()*eigenVector.zy()*eigenVector.zy()*eigenVector.zy();
    const scalar tmp103 = eigenVector.xy()*eigenVector.xy()*eigenVector.xy()*eigenVector.xy();
    const scalar tmp104 = tmp21*tmp61;
    const scalar tmp105 = 6.0*tmp65;
    const scalar tmp106 = tmp26*tmp61;
    const scalar tmp107 = tmp28*tmp63;
    const scalar tmp108 = tmp101*tmp17 + tmp102*tmp18 + tmp103*tmp19 + tmp104*tmp105 + tmp105*tmp107 + 6.0*tmp106*tmp63;
    const scalar tmp109 = tmp21*tmp65;
    const scalar tmp110 = tmp26*tmp63;
    const scalar tmp111 = tmp28*tmp65;
    const scalar tmp112 = tmp63*tmp76;
    const scalar tmp113 = tmp65*tmp78;
    const scalar tmp114 = tmp61*tmp74;
    const scalar tmp115 = 4.0*tmp90;
    const scalar tmp116 = tmp21*tmp86;
    const scalar tmp117 = tmp28*tmp88;
    const scalar tmp118 = tmp104*tmp78 + tmp106*tmp76 + tmp107*tmp78 + tmp109*tmp74 + tmp110*tmp74 + tmp111*tmp76 + tmp112*tmp17 + tmp113*tmp18 + tmp114*tmp19 + tmp115*tmp116 + tmp115*tmp117 + tmp70*tmp86*tmp88;
    const scalar tmp119 = eigenVector.yy()*eigenVector.yy()*eigenVector.yy();
    const scalar tmp120 = tmp119*eigenVector.yx();
    const scalar tmp121 = eigenVector.zy()*eigenVector.zy()*eigenVector.zy();
    const scalar tmp122 = tmp121*eigenVector.zx();
    const scalar tmp123 = eigenVector.xy()*eigenVector.xy()*eigenVector.xy();
    const scalar tmp124 = tmp123*eigenVector.xx();
    const scalar tmp125 = tmp104*tmp42 + tmp106*tmp45 + tmp107*tmp42 + tmp109*tmp40 + tmp110*tmp40 + tmp111*tmp45 + tmp120*tmp17 + tmp122*tmp18 + tmp124*tmp19;
    const scalar tmp126 = 2.0*tmp37;
    const scalar tmp127 = tmp57*tmp63;
    const scalar tmp128 = tmp55*tmp65;
    const scalar tmp129 = tmp53*tmp61;
    const scalar tmp130 = 2.0*tmp90;
    const scalar tmp131 = tmp104*tmp55 + tmp106*tmp57 + tmp107*tmp55 + tmp109*tmp53 + tmp110*tmp53 + tmp111*tmp57 + tmp116*tmp94 + tmp117*tmp94 + tmp127*tmp17 + tmp128*tmp18 + tmp129*tmp19 + tmp130*tmp69 + tmp130*tmp71 + tmp86*tmp97 + tmp88*tmp96;
    const scalar tmp132 = 2.0*tmp52;
    const scalar tmp133 = tmp119*eigenVector.yz();
    const scalar tmp134 = tmp121*eigenVector.zz();
    const scalar tmp135 = tmp123*eigenVector.xz();
    const scalar tmp136 = 3.0*tmp86;
    const scalar tmp137 = 3.0*tmp90;
    const scalar tmp138 = 3.0*tmp88;
    const scalar tmp139 = tmp104*tmp137 + tmp106*tmp138 + tmp107*tmp137 + tmp109*tmp136 + tmp110*tmp136 + tmp111*tmp138 + tmp133*tmp17 + tmp134*tmp18 + tmp135*tmp19;
    const scalar tmp140 = 2.0*tmp92;
    const scalar tmp141 = -tmp20*tmp61*eigenValue.x() - tmp23*tmp65*eigenValue.z() - tmp25*tmp63*eigenValue.y();
    const scalar tmp142 = eigenVector.yz()*eigenVector.yz()*eigenVector.yz()*eigenVector.yz();
    const scalar tmp143 = eigenVector.zz()*eigenVector.zz()*eigenVector.zz()*eigenVector.zz();
    const scalar tmp144 = eigenVector.xz()*eigenVector.xz()*eigenVector.xz()*eigenVector.xz();
    const scalar tmp145 = 6.0*tmp74;
    const scalar tmp146 = tmp21*tmp78;
    const scalar tmp147 = tmp26*tmp76;
    const scalar tmp148 = tmp28*tmp76;
    const scalar tmp149 = tmp142*tmp17 + tmp143*tmp18 + tmp144*tmp19 + tmp145*tmp146 + tmp145*tmp147 + 6.0*tmp148*tmp78;
    const scalar tmp150 = tmp21*tmp74;
    const scalar tmp151 = tmp26*tmp74;
    const scalar tmp152 = tmp28*tmp78;
    const scalar tmp153 = tmp44*tmp76;
    const scalar tmp154 = tmp41*tmp78;
    const scalar tmp155 = tmp39*tmp74;
    const scalar tmp156 = tmp116*tmp93 + tmp117*tmp93 + tmp130*tmp82 + tmp130*tmp83 + tmp146*tmp39 + tmp147*tmp39 + tmp148*tmp41 + tmp150*tmp41 + tmp151*tmp44 + tmp152*tmp44 + tmp153*tmp17 + tmp154*tmp18 + tmp155*tmp19 + tmp53*tmp88*tmp95 + tmp57*tmp86*tmp95;
    const scalar tmp157 = eigenVector.yz()*eigenVector.yz()*eigenVector.yz();
    const scalar tmp158 = tmp157*eigenVector.yx();
    const scalar tmp159 = eigenVector.zz()*eigenVector.zz()*eigenVector.zz();
    const scalar tmp160 = tmp159*eigenVector.zx();
    const scalar tmp161 = eigenVector.xz()*eigenVector.xz()*eigenVector.xz();
    const scalar tmp162 = tmp161*eigenVector.xx();
    const scalar tmp163 = tmp146*tmp54 + tmp147*tmp54 + tmp148*tmp56 + tmp150*tmp56 + tmp151*tmp58 + tmp152*tmp58 + tmp158*tmp17 + tmp160*tmp18 + tmp162*tmp19;
    const scalar tmp164 = tmp157*eigenVector.yy();
    const scalar tmp165 = tmp159*eigenVector.zy();
    const scalar tmp166 = tmp161*eigenVector.xy();
    const scalar tmp167 = tmp136*tmp146 + tmp136*tmp147 + tmp137*tmp148 + tmp137*tmp150 + tmp138*tmp151 + tmp138*tmp152 + tmp164*tmp17 + tmp165*tmp18 + tmp166*tmp19;
    const scalar tmp168 = -tmp20*tmp74*eigenValue.x() - tmp23*tmp78*eigenValue.z() - tmp25*tmp76*eigenValue.y();
    const scalar tmp169 = tmp125*tmp67;
    const scalar tmp170 = tmp156*tmp80;
    const scalar tmp171 = tmp37*tmp72;
    const scalar tmp172 = tmp52*tmp98;
    const scalar tmp173 = tmp131*tmp140;
    const scalar tmp174 = -tmp31*eigenValue.x()*eigenVector.xy() - tmp33*eigenValue.y()*eigenVector.yy() - tmp35*eigenValue.z()*eigenVector.zy();
    const scalar tmp175 = 2.0*D.xy();
    const scalar tmp176 = tmp131*tmp67;
    const scalar tmp177 = tmp163*tmp80;
    const scalar tmp178 = tmp37*tmp98;
    const scalar tmp179 = tmp52*tmp84;
    const scalar tmp180 = tmp140*tmp156;
    const scalar tmp181 = -tmp31*eigenValue.x()*eigenVector.xz() - tmp33*eigenValue.y()*eigenVector.yz() - tmp35*eigenValue.z()*eigenVector.zz();
    const scalar tmp182 = 2.0*D.xz();
    const scalar tmp183 = -tmp20*eigenValue.x()*eigenVector.xy()*eigenVector.xz() - tmp23*eigenValue.z()*eigenVector.zy()*eigenVector.zz() - tmp25*eigenValue.y()*eigenVector.yy()*eigenVector.yz();
    const scalar tmp184 = 2.0*D.yz();
    const scalar tmp185 = tmp120 + tmp122 + tmp124;
    const scalar tmp186 = tmp153 + tmp154 + tmp155;
    const scalar tmp187 = tmp118*tmp186;
    const scalar tmp188 = tmp127 + tmp128 + tmp129;
    const scalar tmp189 = tmp139*tmp188;
    const scalar tmp190 = -tmp119*eigenValue.y()*eigenVector.yx() - tmp121*eigenValue.z()*eigenVector.zx() - tmp123*eigenValue.x()*eigenVector.xx();
    const scalar tmp191 = tmp185*tmp72;
    const scalar tmp192 = tmp186*tmp84;
    const scalar tmp193 = tmp47*tmp67;
    const scalar tmp194 = tmp59*tmp92;
    const scalar tmp195 = 2.0*tmp188;
    const scalar tmp196 = tmp195*tmp98;
    const scalar tmp197 = 2.0*tmp156;
    const scalar tmp198 = -tmp74*eigenValue.x()*eigenVector.xx()*eigenVector.xy() - tmp76*eigenValue.y()*eigenVector.yx()*eigenVector.yy() - tmp78*eigenValue.z()*eigenVector.zx()*eigenVector.zy();
    const scalar tmp199 = tmp125*tmp185;
    const scalar tmp200 = tmp156*tmp186;
    const scalar tmp201 = tmp131*tmp188;
    const scalar tmp202 = 2.0*tmp201;
    const scalar tmp203 = tmp167*tmp186;
    const scalar tmp204 = tmp118*tmp188;
    const scalar tmp205 = -tmp61*eigenValue.x()*eigenVector.xx()*eigenVector.xz() - tmp63*eigenValue.y()*eigenVector.yx()*eigenVector.yz() - tmp65*eigenValue.z()*eigenVector.zx()*eigenVector.zz();
    const scalar tmp206 = tmp131*tmp185;
    const scalar tmp207 = tmp163*tmp186;
    const scalar tmp208 = tmp67*tmp98;
    const scalar tmp209 = tmp84*tmp92;
    const scalar tmp210 = tmp156*tmp195;
    const scalar tmp211 = tmp175*(tmp0*(-tmp100 - tmp141 - tmp199 - tmp200 - tmp202 - tmp48 - 2.0*tmp73) + tmp72) + tmp182*(tmp0*(-tmp183 - tmp206 - tmp207 - 2.0*tmp208 - 2.0*tmp209 - tmp210 - tmp37*tmp59) + tmp98) + tmp184*(tmp0*(-tmp139*tmp185 - 2.0*tmp176 - tmp178 - tmp180 - tmp203 - 2.0*tmp204 - tmp205) + tmp131) + (tmp0*(-tmp174 - tmp191 - tmp192 - 2.0*tmp193 - 2.0*tmp194 - tmp196 - tmp30*tmp37) + tmp47)*D.xx() + (tmp0*(-tmp108*tmp185 - 2.0*tmp169 - tmp171 - tmp173 - tmp187 - 2.0*tmp189 - tmp190) + tmp125)*D.yy() + (tmp0*(-tmp118*tmp185 - tmp140*tmp163 - tmp149*tmp186 - tmp167*tmp195 - tmp197*tmp67 - tmp198 - tmp37*tmp84) + tmp156)*D.zz();
    const scalar tmp212 = tmp158 + tmp160 + tmp162;
    const scalar tmp213 = -tmp157*eigenValue.y()*eigenVector.yx() - tmp159*eigenValue.z()*eigenVector.zx() - tmp161*eigenValue.x()*eigenVector.xx();
    const scalar tmp214 = tmp212*tmp84;
    const scalar tmp215 = tmp188*tmp72;
    const scalar tmp216 = tmp59*tmp80;
    const scalar tmp217 = tmp47*tmp92;
    const scalar tmp218 = 2.0*tmp186;
    const scalar tmp219 = tmp218*tmp98;
    const scalar tmp220 = 2.0*tmp131;
    const scalar tmp221 = tmp163*tmp212;
    const scalar tmp222 = 2.0*tmp200;
    const scalar tmp223 = tmp156*tmp212;
    const scalar tmp224 = tmp125*tmp188;
    const scalar tmp225 = tmp80*tmp98;
    const scalar tmp226 = tmp72*tmp92;
    const scalar tmp227 = tmp131*tmp218;
    const scalar tmp228 = tmp175*(tmp0*(-tmp183 - tmp223 - tmp224 - 2.0*tmp225 - 2.0*tmp226 - tmp227 - tmp47*tmp52) + tmp98) + tmp182*(tmp0*(-tmp100 - tmp168 - tmp201 - tmp221 - tmp222 - tmp60 - 2.0*tmp85) + tmp84) + tmp184*(tmp0*(-tmp167*tmp212 - 2.0*tmp170 - tmp172 - tmp173 - 2.0*tmp187 - tmp189 - tmp198) + tmp156) + (tmp0*(-tmp181 - tmp214 - tmp215 - 2.0*tmp216 - 2.0*tmp217 - tmp219 - tmp30*tmp52) + tmp59)*D.xx() + (tmp0*(-tmp108*tmp188 - tmp118*tmp212 - tmp125*tmp140 - tmp139*tmp218 - tmp205 - tmp220*tmp80 - tmp52*tmp72) + tmp131)*D.yy() + (tmp0*(-tmp149*tmp212 - 2.0*tmp177 - tmp179 - tmp180 - 2.0*tmp203 - tmp204 - tmp213) + tmp163)*D.zz();
    const scalar tmp229 = tmp101 + tmp102 + tmp103;
    const scalar tmp230 = tmp133 + tmp134 + tmp135;
    const scalar tmp231 = tmp139*tmp230;
    const scalar tmp232 = tmp112 + tmp113 + tmp114;
    const scalar tmp233 = tmp118*tmp232;
    const scalar tmp234 = 2.0*tmp185;
    const scalar tmp235 = 2.0*tmp230;
    const scalar tmp236 = -tmp61*tmp74*eigenValue.x() - tmp63*tmp76*eigenValue.y() - tmp65*tmp78*eigenValue.z();
    const scalar tmp237 = tmp156*tmp232;
    const scalar tmp238 = tmp131*tmp230;
    const scalar tmp239 = tmp167*tmp232;
    const scalar tmp240 = tmp118*tmp230;
    const scalar tmp241 = -tmp119*eigenValue.y()*eigenVector.yz() - tmp121*eigenValue.z()*eigenVector.zz() - tmp123*eigenValue.x()*eigenVector.xz();
    const scalar tmp242 = tmp164 + tmp165 + tmp166;
    const scalar tmp243 = -tmp157*eigenValue.y()*eigenVector.yy() - tmp159*eigenValue.z()*eigenVector.zy() - tmp161*eigenValue.x()*eigenVector.xy();
    const scalar tmp244 = tmp118*tmp242;
    const scalar tmp245 = tmp139*tmp232;
    const scalar tmp246 = 2.0*tmp98;
    const scalar tmp247 = tmp167*tmp242;
    const scalar tmp248 = tmp156*tmp242;
    const scalar tmp249 = tmp131*tmp232;
    const scalar tmp250 = tmp175*(tmp0*(-tmp125*tmp230 - tmp205 - 2.0*tmp215 - tmp217 - tmp219 - tmp248 - 2.0*tmp249) + tmp131) + tmp182*(tmp0*(-tmp163*tmp242 - 2.0*tmp192 - tmp194 - tmp196 - tmp198 - 2.0*tmp237 - tmp238) + tmp156) + tmp184*(tmp0*(-tmp202 - tmp222 - tmp231 - 2.0*tmp233 - tmp236 - tmp247 - tmp99) + tmp118) + (tmp0*(-tmp183 - tmp195*tmp47 - tmp218*tmp59 - tmp230*tmp72 - tmp232*tmp246 - tmp242*tmp84 - tmp30*tmp92) + tmp98)*D.xx() + (tmp0*(-tmp108*tmp230 - 2.0*tmp224 - tmp226 - tmp227 - tmp241 - tmp244 - 2.0*tmp245) + tmp139)*D.yy() + (tmp0*(-tmp149*tmp242 - 2.0*tmp207 - tmp209 - tmp210 - 2.0*tmp239 - tmp240 - tmp243) + tmp167)*D.zz();
    const scalar tmp251 = tmp142 + tmp143 + tmp144;


    result.xx() = tmp175*(tmp0*(-tmp169 - tmp170 - 2.0*tmp171 - 2.0*tmp172 - tmp173 - tmp174 - tmp4*tmp47) + tmp47) + tmp182*(tmp0*(-tmp176 - tmp177 - 2.0*tmp178 - 2.0*tmp179 - tmp180 - tmp181 - tmp4*tmp59) + tmp59) + tmp184*(tmp0*(-tmp118*tmp140 - tmp126*tmp131 - tmp132*tmp156 - tmp139*tmp67 - tmp167*tmp80 - tmp183 - tmp4*tmp98) + tmp98) + (tmp0*(-tmp108*tmp67 - tmp118*tmp80 - tmp125*tmp126 - tmp131*tmp132 - tmp139*tmp140 - tmp141 - tmp4*tmp72) + tmp72)*D.yy() + (tmp0*(-tmp118*tmp67 - tmp126*tmp156 - tmp132*tmp163 - tmp140*tmp167 - tmp149*tmp80 - tmp168 - tmp4*tmp84) + tmp84)*D.zz() + (tmp0*(tmp1*eigenValue.x() - tmp100 + tmp2*eigenValue.y() + tmp3*eigenValue.z() - tmp30*tmp4 - 2.0*tmp48 - 2.0*tmp60 - tmp73 - tmp85) + tmp30)*D.xx();
    result.xy() = tmp211;
    result.xz() = tmp228;
    result.yy() = tmp175*(tmp0*(-tmp125*tmp229 - tmp190 - 2.0*tmp191 - tmp193 - tmp196 - tmp237 - 2.0*tmp238) + tmp125) + tmp182*(tmp0*(-tmp131*tmp229 - tmp163*tmp232 - tmp195*tmp84 - tmp197*tmp230 - tmp205 - tmp234*tmp98 - tmp59*tmp67) + tmp131) + tmp184*(tmp0*(-tmp139*tmp229 - 2.0*tmp206 - tmp208 - tmp210 - tmp239 - 2.0*tmp240 - tmp241) + tmp139) + (tmp0*(-tmp141 - tmp195*tmp59 - tmp229*tmp72 - tmp232*tmp84 - tmp234*tmp47 - tmp235*tmp98 - tmp30*tmp67) + tmp72)*D.xx() + (tmp0*(-tmp118*tmp229 - tmp149*tmp232 - tmp163*tmp195 - tmp167*tmp235 - tmp185*tmp197 - tmp236 - tmp67*tmp84) + tmp118)*D.zz() + (tmp0*(tmp101*eigenValue.y() + tmp102*eigenValue.z() + tmp103*eigenValue.x() - tmp108*tmp229 - 2.0*tmp199 - tmp202 - 2.0*tmp231 - tmp233 - tmp73) + tmp108)*D.yy();
    result.yz() = tmp250;
    result.zz() = tmp175*(tmp0*(-tmp125*tmp232 - tmp156*tmp251 - tmp198 - tmp212*tmp246 - tmp218*tmp72 - tmp220*tmp242 - tmp47*tmp80) + tmp156) + tmp182*(tmp0*(-tmp163*tmp251 - tmp213 - 2.0*tmp214 - tmp216 - tmp219 - 2.0*tmp248 - tmp249) + tmp163) + tmp184*(tmp0*(-tmp167*tmp251 - 2.0*tmp223 - tmp225 - tmp227 - tmp243 - 2.0*tmp244 - tmp245) + tmp167) + (tmp0*(-tmp168 - 2.0*tmp212*tmp59 - tmp218*tmp47 - tmp232*tmp72 - tmp242*tmp246 - tmp251*tmp84 - tmp30*tmp80) + tmp84)*D.xx() + (tmp0*(-tmp108*tmp232 - tmp118*tmp251 - tmp125*tmp218 - 2.0*tmp139*tmp242 - tmp212*tmp220 - tmp236 - tmp72*tmp80) + tmp118)*D.yy() + (tmp0*(tmp142*eigenValue.y() + tmp143*eigenValue.z() + tmp144*eigenValue.x() - tmp149*tmp251 - 2.0*tmp221 - tmp222 - tmp233 - 2.0*tmp247 - tmp85) + tmp149)*D.zz();
}
