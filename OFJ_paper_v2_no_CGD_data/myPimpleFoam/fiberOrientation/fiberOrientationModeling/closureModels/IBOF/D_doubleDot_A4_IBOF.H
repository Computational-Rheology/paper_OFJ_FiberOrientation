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
    D. H. Chung and T. H. Kwon, 
    “Invariant-based optimal fitting closure approximation for the numerical prediction of flow-induced fiber orientation,” 
    J. Rheol. (N. Y. N. Y)., vol. 46, no. 1, pp. 169–194, 2002, doi: 10.1122/1.1423312.
    
\*-------------------------------------------------------------------------------------------------------*/

void Foam::fiberOrientation::closureModels::IBOF::D_doubleDot_A4_IBOF
(
    symmTensor& result, 
    const Foam::symmTensor& D, 
    const Foam::symmTensor& A
)
{
    const scalar tmp0 = A.xx()*A.xx();
    const scalar tmp1 = A.xy()*A.xy();
    const scalar tmp2 = A.xz()*A.xz();
    const scalar tmp3 = tmp1 + tmp2;
    const scalar tmp4 = tmp0 + tmp3;
    const scalar tmp5 = 0.5*A.xy();
    const scalar tmp6 = A.xx()*A.xy();
    const scalar tmp7 = A.xy()*A.yy();
    const scalar tmp8 = A.xz()*A.yz();
    const scalar tmp9 = tmp6 + tmp7 + tmp8;
    const scalar tmp10 = 0.5*A.xx();
    const scalar tmp11 = A.yy()*A.yy();
    const scalar tmp12 = A.yz()*A.yz();
    const scalar tmp13 = A.zz()*A.zz();
    const scalar tmp14 = A.xx()*A.yy();
    const scalar tmp15 = tmp14*A.zz();
    const scalar tmp16 = tmp8*A.xy();
    const scalar tmp17 = pow(A.xx() + A.yy() + A.zz(), 2);
    const scalar tmp18 = 0.5*tmp0 + 0.5*tmp11 + tmp12 + 0.5*tmp13 - 0.5*tmp17 + tmp3;
    const scalar tmp19 = -tmp18;
    const scalar tmp20 = tmp19*tmp19;
    const scalar tmp21 = tmp19*tmp19*tmp19;
    const scalar tmp22 = tmp19*tmp19*tmp19*tmp19;
    const scalar tmp23 = tmp19*tmp19*tmp19*tmp19*tmp19;
    const scalar tmp24 = tmp12*A.xx();
    const scalar tmp25 = tmp1*A.zz();
    const scalar tmp26 = tmp2*A.yy();
    const scalar tmp27 = -tmp24 - tmp25 - tmp26 + A.xx()*A.yy()*A.zz() + 2.0*A.xy()*A.xz()*A.yz();
    const scalar tmp28 = tmp27*tmp27;
    const scalar tmp29 = tmp27*tmp27*tmp27;
    const scalar tmp30 = tmp27*tmp27*tmp27*tmp27;
    const scalar tmp31 = tmp27*tmp27*tmp27*tmp27*tmp27;
    const scalar tmp32 = tmp27*(66965.964947122498*tmp0 + 133931.929894245*tmp1 + 66965.964947122498*tmp11 + 133931.929894245*tmp12 + 66965.964947122498*tmp13 - 66965.964947122498*tmp17 + 133931.929894245*tmp2);
    const scalar tmp33 = 1068.7762439282301*tmp0 + 2137.5524878564602*tmp1 + 1068.7762439282301*tmp11 + 2137.5524878564602*tmp12 + 1068.7762439282301*tmp13 - 1068.7762439282301*tmp17 + 2137.5524878564602*tmp2;
    const scalar tmp34 = 1.4*tmp27*tmp33;
    const scalar tmp35 = tmp28*(4958153.453709905*tmp0 + 9916306.9074198101*tmp1 + 4958153.453709905*tmp11 + 9916306.9074198101*tmp12 + 4958153.453709905*tmp13 - 4958153.453709905*tmp17 + 9916306.9074198101*tmp2);
    const scalar tmp36 = 929746.52961154003*tmp0 + 1859493.0592230801*tmp1 + 929746.52961154003*tmp11 + 1859493.0592230801*tmp12 + 929746.52961154003*tmp13 - 929746.52961154003*tmp17 + 1859493.0592230801*tmp2;
    const scalar tmp37 = 1.4*tmp28*tmp36;
    const scalar tmp38 = tmp29*(16850541.013691049*tmp0 + 33701082.027382098*tmp1 + 16850541.013691049*tmp11 + 33701082.027382098*tmp12 + 16850541.013691049*tmp13 - 16850541.013691049*tmp17 + 33701082.027382098*tmp2);
    const scalar tmp39 = 7395951.3822101001*tmp0 + 14791902.7644202*tmp1 + 7395951.3822101001*tmp11 + 14791902.7644202*tmp12 + 7395951.3822101001*tmp13 - 7395951.3822101001*tmp17 + 14791902.7644202*tmp2;
    const scalar tmp40 = 1.4*tmp29*tmp39;
    const scalar tmp41 = -243546726.4462975*tmp0 - 487093452.89259499*tmp1 - 243546726.4462975*tmp11 - 487093452.89259499*tmp12 - 243546726.4462975*tmp13 + 243546726.4462975*tmp17 - 487093452.89259499*tmp2;
    const scalar tmp42 = 1.4*tmp30*tmp41;
    const scalar tmp43 = tmp30*(-785383511.86101997*tmp0 - 1570767023.7220399*tmp1 - 785383511.86101997*tmp11 - 1570767023.7220399*tmp12 - 785383511.86101997*tmp13 + 785383511.86101997*tmp17 - 1570767023.7220399*tmp2);
    const scalar tmp44 = tmp20*tmp27;
    const scalar tmp45 = tmp22*tmp27;
    const scalar tmp46 = tmp20*tmp28;
    const scalar tmp47 = tmp27*(30250.3056757796*tmp0 + 60500.611351559201*tmp1 + 30250.3056757796*tmp11 + 60500.611351559201*tmp12 + 30250.3056757796*tmp13 - 30250.3056757796*tmp17 + 60500.611351559201*tmp2);
    const scalar tmp48 = tmp28*(-2995332.4334491799*tmp0 - 5990664.8668983597*tmp1 - 2995332.4334491799*tmp11 - 5990664.8668983597*tmp12 - 2995332.4334491799*tmp13 + 2995332.4334491799*tmp17 - 5990664.8668983597*tmp2);
    const scalar tmp49 = tmp29*(-283712455.50391853*tmp0 - 567424911.00783706*tmp1 - 283712455.50391853*tmp11 - 567424911.00783706*tmp12 - 283712455.50391853*tmp13 + 283712455.50391853*tmp17 - 567424911.00783706*tmp2);
    const scalar tmp50 = tmp30*(-3962634249.4110899*tmp0 - 7925268498.8221798*tmp1 - 3962634249.4110899*tmp11 - 7925268498.8221798*tmp12 - 3962634249.4110899*tmp13 + 3962634249.4110899*tmp17 - 7925268498.8221798*tmp2);
    const scalar tmp51 = tmp21*tmp27;
    const scalar tmp52 = tmp21*tmp28;
    const scalar tmp53 = tmp20*tmp29;
    const scalar tmp54 = -157.59180543833327*tmp0 + 5843.1014121627168*tmp1*A.zz() - 315.18361087666653*tmp1 - 157.59180543833327*tmp11 + 5843.1014121627168*tmp12*A.xx() - 315.18361087666653*tmp12 - 157.59180543833327*tmp13 - 5843.1014121627168*tmp15 - 11686.202824325434*tmp16 + 157.59180543833327*tmp17 + 5843.1014121627168*tmp2*A.yy() - 315.18361087666653*tmp2 + 6901136.4625648148*tmp20*tmp29 - 2431.652733354043*tmp20 + 1754664.0928393295*tmp21*tmp27 + 48583772.653129004*tmp21*tmp28 + 8607.3235851231639*tmp21 - 11204.663193945982*tmp22 - 760.13765306859932*tmp23 - 873078.28093423136*tmp28 - 9875809.5555594936*tmp29 + 294728034.84386063*tmp30 + 5408425687.0350599*tmp31 - 0.8*tmp32 - tmp34 - 0.8*tmp35 - tmp37 - 0.8*tmp38 - tmp40 - tmp42 - 0.8*tmp43 - 698944.86700303352*tmp44 - 1368943.5740614603*tmp45 - 40573628.86741145*tmp46 - (0.66666666666666667*tmp0 + 1.3333333333333333*tmp1 + 0.66666666666666667*tmp11 + 1.3333333333333333*tmp12 + 0.66666666666666667*tmp13 - 0.66666666666666667*tmp17 + 1.3333333333333333*tmp2 + 1.0)*(247.22882602352038*tmp0 + 494.45765204704077*tmp1 + 247.22882602352038*tmp11 + 494.45765204704077*tmp12 + 247.22882602352038*tmp13 + 6879.1151319721794*tmp15 + 13758.230263944359*tmp16 - 247.22882602352038*tmp17 + 494.45765204704077*tmp2 + 3834.6384047050678*tmp20 - 13278.832221188281*tmp21 + 15476.047042344479*tmp22 + 5601.210975515819*tmp23 - 6879.1151319721794*tmp24 - 6879.1151319721794*tmp25 - 6879.1151319721794*tmp26 - 58225.536367777553*tmp28 - 55265229.681683518*tmp29 - 1833034259.478168*tmp30 - 15366093393.53508*tmp31 - 57260.848802108041*tmp44 - 5991860.9531104071*tmp45 - 66792738.808180198*tmp46 + 1.2*tmp47 + 1.2*tmp48 + 1.2*tmp49 + 1.2*tmp50 + 2436515.5238744877*tmp51 + 158549793.7719996*tmp52 - 1948319935.4517961*tmp53 + 28.097554988519878) - 19.256621626473745;
    const scalar tmp55 = 206.02402168626699*tmp0 + 412.04804337253398*tmp1 + 206.02402168626699*tmp11 + 412.04804337253398*tmp12 + 206.02402168626699*tmp13 + 5732.5959433101498*tmp15 + 11465.1918866203*tmp16 - 206.02402168626699*tmp17 + 412.04804337253398*tmp2 + 3195.53200392089*tmp20 - 11065.6935176569*tmp21 + 12896.7058686204*tmp22 + 4667.6758129298496*tmp23 - 5732.5959433101498*tmp24 - 5732.5959433101498*tmp25 - 5732.5959433101498*tmp26 - 48521.280306481298*tmp28 - 46054358.0680696*tmp29 - 1527528549.56514*tmp30 - 12805077827.9459*tmp31 - 47717.374001756703*tmp44 - 4993217.4609253397*tmp45 - 55660615.673483498*tmp46 + tmp47 + tmp48 + tmp49 + tmp50 + 2030429.6032287399*tmp51 + 132124828.143333*tmp52 - 1623599946.20983*tmp53 + 23.414629157099899;
    const scalar tmp56 = tmp4*tmp55;
    const scalar tmp57 = 217.55057658016449*tmp0 + 435.10115316032898*tmp1 + 217.55057658016449*tmp11 + 435.10115316032898*tmp12 + 217.55057658016449*tmp13 + 7034.4365791647597*tmp15 + 14068.873158329519*tmp16 - 217.55057658016449*tmp17 + 435.10115316032898*tmp2 + 3723.89335663877*tmp20 - 15939.2396237307*tmp21 + 32221.941625641699*tmp22 - 23215.348852529802*tmp23 - 7034.4365791647597*tmp24 - 7034.4365791647597*tmp25 - 7034.4365791647597*tmp26 + 823995.18736610597*tmp28 + 8009700.2684979597*tmp29 - 257258805.87056699*tmp30 - 3957693983.0447302*tmp31 + tmp32 + tmp35 + tmp38 + tmp43 + 880683.51532791601*tmp44 + 2144190.90344474*tmp45 + 37901059.935526699*tmp46 - 2370104.5868925201*tmp51 - 44927559.185148999*tmp52 - 21313392.0223355*tmp53 + 24.940908165785999;
    const scalar tmp58 = -11.749039875570251*tmp0 - 23.498079751140502*tmp1 - 11.749039875570251*tmp11 - 23.498079751140502*tmp12 - 11.749039875570251*tmp13 + 11.749039875570251*tmp17 - 23.498079751140502*tmp2 + 7323414.9421357801*tmp20*tmp28 + 7249697.9680739902*tmp20*tmp29 - 391.04425139783802*tmp20 + 101013.983339062*tmp21*tmp27 + 2960.0486527581402*tmp21 - 10409.2072189767*tmp22 + 13808.8690964946*tmp23 - 153.96582059350601*tmp24 - 153.96582059350601*tmp25 - 153.96582059350601*tmp26 + tmp27*tmp33 + tmp28*tmp36 + 152772.950743819*tmp28 + tmp29*tmp39 + 2477178.1005436601*tmp29 + tmp30*tmp41 - 63514992.962433599*tmp30 - 1601621786.1423399*tmp31 - 4001.3894709281199*tmp44 - 247435.10621023699*tmp45 - 9029803.7892927192*tmp52 + 153.96582059350601*A.xx()*A.yy()*A.zz() + 307.93164118701202*A.xy()*A.xz()*A.yz() - 0.49721779011075401;
    const scalar tmp59 = 43.510115316032902*tmp0 + 87.020230632065804*tmp1 + 43.510115316032902*tmp11 + 87.020230632065804*tmp12 + 43.510115316032902*tmp13 + 1406.887315832952*tmp15 + 2813.7746316659041*tmp16 - 43.510115316032902*tmp17 + 87.020230632065804*tmp2 + 744.77867132775407*tmp20 - 3187.8479247461401*tmp21 + 6444.3883251283405*tmp22 - 4643.0697705059602*tmp23 - 1406.887315832952*tmp24 - 1406.887315832952*tmp25 - 1406.887315832952*tmp26 + 164799.03747322119*tmp28 + 1601940.0536995921*tmp29 - 51451761.1741134*tmp30 - 791538796.60894608*tmp31 + 0.2*tmp32 + 0.2*tmp35 + 0.2*tmp38 + 0.2*tmp43 + 176136.70306558322*tmp44 + 428838.18068894802*tmp45 + 7580211.9871053398*tmp46 - 474020.91737850406*tmp51 - 8985511.8370297998*tmp52 - 4262678.4044671003*tmp53 + 4.9881816331572004;
    const scalar tmp60 = -0.85714285714285714*tmp55*(-0.4*tmp0 - 0.8*tmp1 - 0.4*tmp11 - 0.8*tmp12 - 0.4*tmp13 + 0.4*tmp17 - 0.8*tmp2 - 1.6*tmp20 - 0.66666666666666667*tmp24 - 0.66666666666666667*tmp25 - 0.66666666666666667*tmp26 + 0.66666666666666667*A.xx()*A.yy()*A.zz() + 1.3333333333333333*A.xy()*A.xz()*A.yz() - 0.2) - 0.85714285714285714*tmp59*(-2.0*tmp0 - 4.0*tmp1 - 2.0*tmp11 - 4.0*tmp12 - 2.0*tmp13 + 2.0*tmp17 - 4.0*tmp2 + 1.0) + 0.85714285714285714*(tmp18 + 0.16666666666666667)*(-16.448655825798351*tmp0 - 32.897311651596702*tmp1 - 16.448655825798351*tmp11 - 32.897311651596702*tmp12 - 16.448655825798351*tmp13 + 16.448655825798351*tmp17 - 32.897311651596702*tmp2 + 10252780.918990092*tmp20*tmp28 + 10149577.155303586*tmp20*tmp29 - 547.46195195697317*tmp20 + 141419.5766746868*tmp21*tmp27 + 4144.0681138613963*tmp21 - 14572.89010656738*tmp22 + 19332.41673509244*tmp23 - 215.55214883090841*tmp24 - 215.55214883090841*tmp25 - 215.55214883090841*tmp26 + 213882.13104134658*tmp28 + 3468049.3407611242*tmp29 - 88920990.14740704*tmp30 - 2242270500.5992756*tmp31 + tmp34 + tmp37 + tmp40 + tmp42 - 5601.9452592993675*tmp44 - 346409.14869433176*tmp45 - 12641725.305009807*tmp52 + 215.55214883090841*A.xx()*A.yy()*A.zz() + 431.10429766181682*A.xy()*A.xz()*A.yz() - 0.6961049061550556) + 0.85714285714285714;
    const scalar tmp61 = tmp5*tmp60 + tmp58*(0.5*tmp6 + 0.5*tmp7 + 0.5*tmp8);
    const scalar tmp62 = tmp54*(tmp10*tmp9 + tmp4*tmp5) + tmp56*tmp9 + tmp57*tmp6 + tmp61;
    const scalar tmp63 = 2.0*D.xy();
    const scalar tmp64 = 0.5*A.xz();
    const scalar tmp65 = A.xx()*A.xz();
    const scalar tmp66 = A.xy()*A.yz();
    const scalar tmp67 = A.xz()*A.zz();
    const scalar tmp68 = tmp65 + tmp66 + tmp67;
    const scalar tmp69 = tmp58*(0.5*tmp65 + 0.5*tmp66 + 0.5*tmp67) + tmp60*tmp64;
    const scalar tmp70 = tmp54*(tmp10*tmp68 + tmp4*tmp64) + tmp56*tmp68 + tmp57*tmp65 + tmp69;
    const scalar tmp71 = 2.0*D.xz();
    const scalar tmp72 = 0.33333333333333333*A.xx();
    const scalar tmp73 = A.xy()*A.xz();
    const scalar tmp74 = 0.16666666666666667*A.yy();
    const scalar tmp75 = 0.16666666666666667*A.zz();
    const scalar tmp76 = A.yy()*A.yz();
    const scalar tmp77 = A.yz()*A.zz();
    const scalar tmp78 = tmp73 + tmp76 + tmp77;
    const scalar tmp79 = 0.33333333333333333*tmp4;
    const scalar tmp80 = 0.66666666666666667*tmp9;
    const scalar tmp81 = 0.16666666666666667*A.yz();
    const scalar tmp82 = 0.33333333333333333*tmp9;
    const scalar tmp83 = 0.33333333333333333*A.xy();
    const scalar tmp84 = 0.16666666666666667*A.xx();
    const scalar tmp85 = tmp54*(tmp4*tmp81 + tmp68*tmp83 + tmp78*tmp84 + tmp82*A.xz()) + tmp55*(tmp68*tmp80 + tmp78*tmp79) + tmp57*(tmp72*A.yz() + 0.66666666666666667*tmp73) + tmp58*(0.16666666666666667*tmp73 + tmp74*A.yz() + tmp75*A.yz()) + tmp60*tmp81;
    const scalar tmp86 = 2.0*D.yz();
    const scalar tmp87 = tmp58*(0.26666666666666667*tmp0 + 0.53333333333333333*tmp1 + 0.26666666666666667*tmp11 + 0.53333333333333333*tmp12 + 0.26666666666666667*tmp13 - 0.93333333333333333*tmp15 - 1.8666666666666667*tmp16 - 0.26666666666666667*tmp17 + 0.53333333333333333*tmp2 + 0.93333333333333333*tmp24 + 0.93333333333333333*tmp25 + 0.93333333333333333*tmp26 + 0.2);
    const scalar tmp88 = tmp55*(0.057142857142857143*tmp0 + 0.11428571428571429*tmp1 + 0.057142857142857143*tmp11 + 0.11428571428571429*tmp12 + 0.057142857142857143*tmp13 - 0.22857142857142857*tmp15 - 0.45714285714285714*tmp16 - 0.057142857142857143*tmp17 + 0.11428571428571429*tmp2 + 0.22857142857142857*tmp20 + 0.22857142857142857*tmp24 + 0.22857142857142857*tmp25 + 0.22857142857142857*tmp26 + tmp27*(-0.53333333333333333*tmp0 - 1.0666666666666667*tmp1 - 0.53333333333333333*tmp11 - 1.0666666666666667*tmp12 - 0.53333333333333333*tmp13 + 0.53333333333333333*tmp17 - 1.0666666666666667*tmp2) + 0.028571428571428571);
    const scalar tmp89 = tmp59*(-0.28571428571428571*tmp0 - 0.57142857142857143*tmp1 - 0.28571428571428571*tmp11 - 0.57142857142857143*tmp12 - 0.28571428571428571*tmp13 + 0.28571428571428571*tmp17 - 0.57142857142857143*tmp2 - 2.6666666666666667*tmp24 - 2.6666666666666667*tmp25 - 2.6666666666666667*tmp26 + 2.6666666666666667*A.xx()*A.yy()*A.zz() + 5.3333333333333333*A.xy()*A.xz()*A.yz() + 0.14285714285714286);
    const scalar tmp90 = -0.6*tmp87 - 0.6*tmp88 + 0.6*tmp89 - 0.085714285714285714;
    const scalar tmp91 = tmp1 + tmp11 + tmp12;
    const scalar tmp92 = 0.16666666666666667*tmp0 + 0.16666666666666667*tmp12;
    const scalar tmp93 = 0.16666666666666667*tmp11 + 0.16666666666666667*tmp2;
    const scalar tmp94 = -0.2*tmp87 - 0.2*tmp88 + 0.2*tmp89 - 0.028571428571428571;
    const scalar tmp95 = tmp54*(tmp4*tmp74 + tmp80*A.xy() + tmp84*tmp91) + tmp55*(tmp79*tmp91 + 0.66666666666666667*(tmp9*tmp9)) + tmp57*(0.66666666666666667*tmp1 + 0.33333333333333333*tmp14) + tmp58*(0.33333333333333333*tmp1 + tmp92 + tmp93) + tmp60*(tmp74 + tmp84) + tmp94;
    const scalar tmp96 = tmp12 + tmp13 + tmp2;
    const scalar tmp97 = 0.66666666666666667*tmp68;
    const scalar tmp98 = 0.16666666666666667*tmp1 + 0.16666666666666667*tmp13;
    const scalar tmp99 = tmp54*(tmp4*tmp75 + tmp84*tmp96 + tmp97*A.xz()) + tmp55*(tmp79*tmp96 + 0.66666666666666667*(tmp68*tmp68)) + tmp57*(0.66666666666666667*tmp2 + tmp72*A.zz()) + tmp58*(0.33333333333333333*tmp2 + tmp92 + tmp98) + tmp60*(tmp75 + tmp84) + tmp94;
    const scalar tmp100 = 0.5*A.yy();
    const scalar tmp101 = tmp55*tmp91;
    const scalar tmp102 = tmp101*tmp9 + tmp54*(tmp100*tmp9 + tmp5*tmp91) + tmp57*tmp7 + tmp61;
    const scalar tmp103 = 0.16666666666666667*A.xy();
    const scalar tmp104 = 0.33333333333333333*tmp68;
    const scalar tmp105 = 0.33333333333333333*A.xz();
    const scalar tmp106 = tmp103*tmp60 + tmp54*(tmp103*tmp96 + tmp104*A.yz() + tmp105*tmp78 + tmp75*tmp9) + tmp55*(tmp78*tmp97 + tmp82*tmp96) + tmp57*(0.66666666666666667*tmp8 + tmp83*A.zz()) + tmp58*(tmp74*A.xy() + 0.16666666666666667*tmp8 + tmp84*A.xy());
    const scalar tmp107 = 0.16666666666666667*A.xz();
    const scalar tmp108 = tmp107*tmp60 + tmp54*(tmp107*tmp91 + tmp68*tmp74 + tmp78*tmp83 + tmp82*A.yz()) + tmp55*(tmp104*tmp91 + tmp78*tmp80) + tmp57*(tmp105*A.yy() + 0.66666666666666667*tmp66) + tmp58*(0.16666666666666667*tmp66 + tmp75*A.xz() + tmp84*A.xz());
    const scalar tmp109 = tmp102*D.yy() + tmp106*D.zz() + tmp108*tmp86 + tmp62*D.xx() + tmp63*tmp95 + tmp71*tmp85;
    const scalar tmp110 = 0.5*A.zz();
    const scalar tmp111 = tmp55*tmp96;
    const scalar tmp112 = tmp111*tmp68 + tmp54*(tmp110*tmp68 + tmp64*tmp96) + tmp57*tmp67 + tmp69;
    const scalar tmp113 = tmp106*tmp86 + tmp108*D.yy() + tmp112*D.zz() + tmp63*tmp85 + tmp70*D.xx() + tmp71*tmp99;
    const scalar tmp114 = 0.5*A.yz();
    const scalar tmp115 = tmp114*tmp60 + tmp58*(0.5*tmp73 + 0.5*tmp76 + 0.5*tmp77);
    const scalar tmp116 = tmp101*tmp78 + tmp115 + tmp54*(tmp100*tmp78 + tmp114*tmp91) + tmp57*tmp76;
    const scalar tmp117 = tmp54*(tmp74*tmp96 + tmp75*tmp91 + 0.66666666666666667*tmp78*A.yz()) + tmp55*(0.33333333333333333*tmp91*tmp96 + 0.66666666666666667*(tmp78*tmp78)) + tmp57*(0.66666666666666667*tmp12 + 0.33333333333333333*A.yy()*A.zz()) + tmp58*(0.33333333333333333*tmp12 + tmp93 + tmp98) + tmp60*(tmp74 + tmp75) + tmp94;
    const scalar tmp118 = tmp111*tmp78 + tmp115 + tmp54*(tmp110*tmp78 + tmp114*tmp96) + tmp57*tmp77;
    const scalar tmp119 = tmp106*tmp71 + tmp108*tmp63 + tmp116*D.yy() + tmp117*tmp86 + tmp118*D.zz() + tmp85*D.xx();

    result.xx() = tmp62*tmp63 + tmp70*tmp71 + tmp85*tmp86 + tmp95*D.yy() + tmp99*D.zz() + D.xx()*(tmp0*tmp57 + tmp4*tmp54*A.xx() + tmp4*tmp58 + tmp55*(tmp4*tmp4) + tmp60*A.xx() + tmp90);
    result.xy() = tmp109;
    result.xz() = tmp113;
    result.yy() = tmp102*tmp63 + tmp108*tmp71 + tmp116*tmp86 + tmp117*D.zz() + tmp95*D.xx() + D.yy()*(tmp11*tmp57 + tmp54*tmp91*A.yy() + tmp55*(tmp91*tmp91) + tmp58*tmp91 + tmp60*A.yy() + tmp90);
    result.yz() = tmp119;
    result.zz() = tmp106*tmp63 + tmp112*tmp71 + tmp117*D.yy() + tmp118*tmp86 + tmp99*D.xx() + D.zz()*(tmp13*tmp57 + tmp54*tmp96*A.zz() + tmp55*(tmp96*tmp96) + tmp58*tmp96 + tmp60*A.zz() + tmp90);
}
