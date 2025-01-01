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


void Foam::fiberOrientation::closureModels::ORE::D_doubleDot_A4_ORE
(
    symmTensor& result, 
    const Foam::symmTensor& D,
    const Foam::vector& eigenValue,
    const Foam::tensor& eigenVector
)
{
    const scalar tmp0 = eigenValue.y()*eigenValue.y();
    const scalar tmp1 = eigenValue.y()*eigenValue.y()*eigenValue.y();
    const scalar tmp2 = eigenValue.y()*eigenValue.y()*eigenValue.y()*eigenValue.y();
    const scalar tmp3 = eigenValue.z()*eigenValue.z();
    const scalar tmp4 = eigenValue.z()*eigenValue.z()*eigenValue.z();
    const scalar tmp5 = eigenValue.z()*eigenValue.z()*eigenValue.z()*eigenValue.z();
    const scalar tmp6 = eigenValue.y()*eigenValue.z();
    const scalar tmp7 = tmp3*eigenValue.y();
    const scalar tmp8 = tmp4*eigenValue.y();
    const scalar tmp9 = tmp0*eigenValue.z();
    const scalar tmp10 = tmp1*eigenValue.z();
    const scalar tmp11 = tmp0*tmp3;
    const scalar tmp12 = 8.4367774677832497*tmp0 - 8.6389141928401596*tmp1 + 7.4814687062444101*tmp10 + 9.3252034345266104*tmp11 + 3.59772251134254*tmp2 + 6.8815395205804402*tmp3 - 6.4872893364192601*tmp4 + 2.28476531637958*tmp5 + 11.8273285968852*tmp6 - 15.9120667157641*tmp7 + 7.74683751713295*tmp8 - 15.1515872606307*tmp9 - 3.0370993982540599*eigenValue.y() - 3.3152722974214601*eigenValue.z() + 0.63625679688068704;
    const scalar tmp13 = 11.342092427815899*tmp0 - 12.3875632855619*tmp1 + 11.749291117702599*tmp10 + 9.8159838971674809*tmp11 + 4.8836659777148901*tmp2 + 3.84459692420086*tmp3 - 2.1162321447100401*tmp4 + 0.50804138736663695*tmp5 + 11.958956233232*tmp6 - 10.9582626069691*tmp7 + 3.4790151056743901*tmp8 - 20.727799468413199*tmp9 - 4.4797087319373796*eigenValue.y() - 1.87266296373814*eigenValue.z() + 0.63625679688068704;
    const scalar tmp14 = 25.868475525388401*tmp0 - 26.963691523971601*tmp1 + 26.113491400537502*tmp10 + 27.3346798054488*tmp11 + 10.611741806606*tmp2 + 13.829469912194*tmp3 - 10.8801761133174*tmp4 + 3.43213840334779*tmp5 + 34.319901891698699*tmp6 - 37.702911802938402*tmp7 + 15.2650686148651*tmp8 - 50.275643192748497*tmp9 - 12.2570587036254*eigenValue.y() - 9.1219650978269193*eigenValue.z() + 2.7405328956025299;
    const scalar tmp15 = eigenVector.xx()*eigenVector.xx();
    const scalar tmp16 = 27.925927700265497*tmp0*eigenValue.z() - 14.386895242710526*tmp0 + 15.35617030834667*tmp1 - 15.190656905997844*tmp10 - 13.912730134044835*tmp11 - 5.9488426364891751*tmp2 + 16.374553847071702*tmp3*eigenValue.y() - 5.3962636579072099*tmp3 + 3.2545594608040904*tmp4 - 0.82770723716742345*tmp5 - 17.225764764022749*tmp6 - 5.4986231017032701*tmp8 + 5.8498340186543603*eigenValue.y() + 3.8396778820717996*eigenValue.z() - 0.87026644780126494;
    const scalar tmp17 = tmp15*tmp16;
    const scalar tmp18 = eigenVector.zx()*eigenVector.zx();
    const scalar tmp19 = 6.0*tmp18;
    const scalar tmp20 = eigenVector.yx()*eigenVector.yx();
    const scalar tmp21 = 22.349715492483*tmp0*eigenValue.z() - 11.481580282677875*tmp0 + 11.607521215624931*tmp1 - 10.922834494539657*tmp10 - 13.421949671403965*tmp11 - 4.662899170116825*tmp2 + 21.3283579558667*tmp3*eigenValue.y() - 8.4332062542867909*tmp3 + 7.62561665251331*tmp4 - 2.6044311661803663*tmp5 - 17.09413712767595*tmp6 - 9.7664455131618304*tmp8 + 5.4072246849710401*eigenValue.y() + 4.2822872157551197*eigenValue.z() - 0.87026644780126494;
    const scalar tmp22 = tmp15*tmp21;
    const scalar tmp23 = 3.0448028148946253*tmp0 - 2.9686070227847718*tmp1 + 3.4413657882952453*tmp10 + 4.0967462368773546*tmp11 + 1.065176658774285*tmp2 + 1.5516667337063499*tmp3 - 1.1383273160940499*tmp4 + 0.31966584980078649*tmp5 + 5.2668085307907493*tmp6 - 5.416291240102602*tmp7 + 2.0196079960288804*tmp8 - 7.1981282318522979*tmp9 - 1.3701252867169806*eigenValue.y() - 0.9670149183336596*eigenValue.z() + 0.23400965092057779;
    const scalar tmp24 = tmp20*tmp23;
    const scalar tmp25 = tmp12*(eigenVector.yx()*eigenVector.yx()*eigenVector.yx());
    const scalar tmp26 = tmp13*(eigenVector.zx()*eigenVector.zx()*eigenVector.zx());
    const scalar tmp27 = tmp14*(eigenVector.xx()*eigenVector.xx()*eigenVector.xx());
    const scalar tmp28 = tmp16*tmp18;
    const scalar tmp29 = eigenVector.xx()*eigenVector.xy();
    const scalar tmp30 = 3.0*tmp29;
    const scalar tmp31 = eigenVector.zx()*eigenVector.zy();
    const scalar tmp32 = 3.0*tmp31;
    const scalar tmp33 = tmp20*tmp21;
    const scalar tmp34 = eigenVector.yx()*eigenVector.yy();
    const scalar tmp35 = 3.0*tmp34;
    const scalar tmp36 = tmp18*tmp23;
    const scalar tmp37 = tmp17*tmp32 + tmp22*tmp35 + tmp24*tmp32 + tmp25*eigenVector.yy() + tmp26*eigenVector.zy() + tmp27*eigenVector.xy() + tmp28*tmp30 + tmp30*tmp33 + tmp35*tmp36;
    const scalar tmp38 = 2.0*D.xy();
    const scalar tmp39 = eigenVector.xx()*eigenVector.xz();
    const scalar tmp40 = 3.0*tmp39;
    const scalar tmp41 = eigenVector.zx()*eigenVector.zz();
    const scalar tmp42 = 3.0*tmp41;
    const scalar tmp43 = eigenVector.yx()*eigenVector.yz();
    const scalar tmp44 = 3.0*tmp43;
    const scalar tmp45 = tmp17*tmp42 + tmp22*tmp44 + tmp24*tmp42 + tmp25*eigenVector.yz() + tmp26*eigenVector.zz() + tmp27*eigenVector.xz() + tmp28*tmp40 + tmp33*tmp40 + tmp36*tmp44;
    const scalar tmp46 = 2.0*D.xz();
    const scalar tmp47 = eigenVector.zy()*eigenVector.zy();
    const scalar tmp48 = eigenVector.xy()*eigenVector.xy();
    const scalar tmp49 = eigenVector.yy()*eigenVector.yy();
    const scalar tmp50 = tmp12*tmp20;
    const scalar tmp51 = tmp13*tmp18;
    const scalar tmp52 = tmp14*tmp15;
    const scalar tmp53 = 4.0*tmp31;
    const scalar tmp54 = tmp16*tmp29;
    const scalar tmp55 = 4.0*tmp21;
    const scalar tmp56 = tmp23*tmp34;
    const scalar tmp57 = tmp17*tmp47 + tmp22*tmp49 + tmp24*tmp47 + tmp28*tmp48 + tmp29*tmp34*tmp55 + tmp33*tmp48 + tmp36*tmp49 + tmp47*tmp51 + tmp48*tmp52 + tmp49*tmp50 + tmp53*tmp54 + tmp53*tmp56;
    const scalar tmp58 = eigenVector.zz()*eigenVector.zz();
    const scalar tmp59 = eigenVector.xz()*eigenVector.xz();
    const scalar tmp60 = eigenVector.yz()*eigenVector.yz();
    const scalar tmp61 = 4.0*tmp41;
    const scalar tmp62 = tmp16*tmp39;
    const scalar tmp63 = tmp23*tmp43;
    const scalar tmp64 = tmp17*tmp58 + tmp22*tmp60 + tmp24*tmp58 + tmp28*tmp59 + tmp33*tmp59 + tmp36*tmp60 + tmp39*tmp43*tmp55 + tmp50*tmp60 + tmp51*tmp58 + tmp52*tmp59 + tmp61*tmp62 + tmp61*tmp63;
    const scalar tmp65 = eigenVector.zy()*eigenVector.zz();
    const scalar tmp66 = eigenVector.xy()*eigenVector.xz();
    const scalar tmp67 = eigenVector.yy()*eigenVector.yz();
    const scalar tmp68 = 2.0*tmp41;
    const scalar tmp69 = 2.0*tmp31;
    const scalar tmp70 = 2.0*tmp21;
    const scalar tmp71 = tmp43*tmp70;
    const scalar tmp72 = tmp39*tmp70;
    const scalar tmp73 = tmp17*tmp65 + tmp22*tmp67 + tmp24*tmp65 + tmp28*tmp66 + tmp29*tmp71 + tmp33*tmp66 + tmp34*tmp72 + tmp36*tmp67 + tmp50*tmp67 + tmp51*tmp65 + tmp52*tmp66 + tmp54*tmp68 + tmp56*tmp68 + tmp62*tmp69 + tmp63*tmp69;
    const scalar tmp74 = 2.0*D.yz();
    const scalar tmp75 = eigenVector.yy()*eigenVector.yy()*eigenVector.yy();
    const scalar tmp76 = tmp12*eigenVector.yx();
    const scalar tmp77 = eigenVector.zy()*eigenVector.zy()*eigenVector.zy();
    const scalar tmp78 = tmp13*eigenVector.zx();
    const scalar tmp79 = eigenVector.xy()*eigenVector.xy()*eigenVector.xy();
    const scalar tmp80 = tmp14*eigenVector.xx();
    const scalar tmp81 = tmp16*tmp47;
    const scalar tmp82 = tmp16*tmp48;
    const scalar tmp83 = tmp21*tmp49;
    const scalar tmp84 = tmp21*tmp48;
    const scalar tmp85 = tmp23*tmp47;
    const scalar tmp86 = tmp23*tmp49;
    const scalar tmp87 = tmp30*tmp81 + tmp30*tmp83 + tmp32*tmp82 + tmp32*tmp86 + tmp35*tmp84 + tmp35*tmp85 + tmp75*tmp76 + tmp77*tmp78 + tmp79*tmp80;
    const scalar tmp88 = tmp16*tmp58;
    const scalar tmp89 = tmp16*tmp59;
    const scalar tmp90 = tmp21*tmp60;
    const scalar tmp91 = tmp21*tmp59;
    const scalar tmp92 = tmp23*tmp58;
    const scalar tmp93 = tmp23*tmp60;
    const scalar tmp94 = 2.0*tmp65;
    const scalar tmp95 = tmp16*tmp66;
    const scalar tmp96 = tmp23*tmp67;
    const scalar tmp97 = tmp29*tmp88 + tmp29*tmp90 + tmp31*tmp89 + tmp31*tmp93 + tmp34*tmp91 + tmp34*tmp92 + tmp58*tmp78*eigenVector.zy() + tmp59*tmp80*eigenVector.xy() + tmp60*tmp76*eigenVector.yy() + tmp62*tmp94 + tmp63*tmp94 + tmp66*tmp71 + tmp67*tmp72 + tmp68*tmp95 + tmp68*tmp96;
    const scalar tmp98 = tmp29*tmp67*tmp70 + tmp34*tmp66*tmp70 + tmp39*tmp81 + tmp39*tmp83 + tmp41*tmp82 + tmp41*tmp86 + tmp43*tmp84 + tmp43*tmp85 + tmp47*tmp78*eigenVector.zz() + tmp48*tmp80*eigenVector.xz() + tmp49*tmp76*eigenVector.yz() + tmp54*tmp94 + tmp56*tmp94 + tmp69*tmp95 + tmp69*tmp96;
    const scalar tmp99 = tmp37*D.xx() + tmp38*tmp57 + tmp46*tmp73 + tmp74*tmp98 + tmp87*D.yy() + tmp97*D.zz();
    const scalar tmp100 = eigenVector.yz()*eigenVector.yz()*eigenVector.yz();
    const scalar tmp101 = eigenVector.zz()*eigenVector.zz()*eigenVector.zz();
    const scalar tmp102 = eigenVector.xz()*eigenVector.xz()*eigenVector.xz();
    const scalar tmp103 = tmp100*tmp76 + tmp101*tmp78 + tmp102*tmp80 + tmp40*tmp88 + tmp40*tmp90 + tmp42*tmp89 + tmp42*tmp93 + tmp44*tmp91 + tmp44*tmp92;
    const scalar tmp104 = tmp103*D.zz() + tmp38*tmp73 + tmp45*D.xx() + tmp46*tmp64 + tmp74*tmp97 + tmp98*D.yy();
    const scalar tmp105 = 6.0*tmp47;
    const scalar tmp106 = 3.0*tmp66;
    const scalar tmp107 = 3.0*tmp65;
    const scalar tmp108 = 3.0*tmp67;
    const scalar tmp109 = tmp106*tmp81 + tmp106*tmp83 + tmp107*tmp82 + tmp107*tmp86 + tmp108*tmp84 + tmp108*tmp85 + tmp12*tmp75*eigenVector.yz() + tmp13*tmp77*eigenVector.zz() + tmp14*tmp79*eigenVector.xz();
    const scalar tmp110 = 4.0*tmp65;
    const scalar tmp111 = tmp110*tmp95 + tmp110*tmp96 + tmp12*tmp49*tmp60 + tmp13*tmp47*tmp58 + tmp14*tmp48*tmp59 + tmp47*tmp89 + tmp47*tmp93 + tmp48*tmp88 + tmp48*tmp90 + tmp49*tmp91 + tmp49*tmp92 + tmp55*tmp66*tmp67;
    const scalar tmp112 = tmp100*tmp12*eigenVector.yy() + tmp101*tmp13*eigenVector.zy() + tmp102*tmp14*eigenVector.xy() + tmp106*tmp88 + tmp106*tmp90 + tmp107*tmp89 + tmp107*tmp93 + tmp108*tmp91 + tmp108*tmp92;
    const scalar tmp113 = tmp109*D.yy() + tmp111*tmp74 + tmp112*D.zz() + tmp38*tmp98 + tmp46*tmp97 + tmp73*D.xx();
    const scalar tmp114 = 6.0*tmp59;


    result.xx() = tmp37*tmp38 + tmp45*tmp46 + tmp57*D.yy() + tmp64*D.zz() + tmp73*tmp74 + D.xx()*(tmp12*(eigenVector.yx()*eigenVector.yx()*eigenVector.yx()*eigenVector.yx()) + tmp13*(eigenVector.zx()*eigenVector.zx()*eigenVector.zx()*eigenVector.zx()) + tmp14*(eigenVector.xx()*eigenVector.xx()*eigenVector.xx()*eigenVector.xx()) + tmp17*tmp19 + tmp19*tmp24 + 6.0*tmp20*tmp22);
    result.xy() = tmp99;
    result.xz() = tmp104;
    result.yy() = tmp109*tmp74 + tmp111*D.zz() + tmp38*tmp87 + tmp46*tmp98 + tmp57*D.xx() + D.yy()*(tmp105*tmp82 + tmp105*tmp86 + tmp12*(eigenVector.yy()*eigenVector.yy()*eigenVector.yy()*eigenVector.yy()) + tmp13*(eigenVector.zy()*eigenVector.zy()*eigenVector.zy()*eigenVector.zy()) + tmp14*(eigenVector.xy()*eigenVector.xy()*eigenVector.xy()*eigenVector.xy()) + 6.0*tmp49*tmp84);
    result.yz() = tmp113;
    result.zz() = tmp103*tmp46 + tmp111*D.yy() + tmp112*tmp74 + tmp38*tmp97 + tmp64*D.xx() + D.zz()*(tmp114*tmp88 + tmp114*tmp90 + tmp12*(eigenVector.yz()*eigenVector.yz()*eigenVector.yz()*eigenVector.yz()) + tmp13*(eigenVector.zz()*eigenVector.zz()*eigenVector.zz()*eigenVector.zz()) + tmp14*(eigenVector.xz()*eigenVector.xz()*eigenVector.xz()*eigenVector.xz()) + 6.0*tmp60*tmp92);
}
