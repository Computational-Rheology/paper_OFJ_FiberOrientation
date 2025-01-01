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

   J. H. Phelps and C. L. Tucker, 
   “An anisotropic rotary diffusion model for fiber orientation in short- and long-fiber thermoplastics,”
    J. Nonnewton. Fluid Mech., vol. 156, no. 3, pp. 165–176, 2009, doi: 10.1016/j.jnnfm.2008.08.002.

\*-------------------------------------------------------------------------------------------------------*/

void Foam::fiberOrientation::closureModels::IBOF::RSC_D_doubleDot_A4_IBOF
(
    Foam::symmTensor& result,
    const Foam::symmTensor& D,
    const Foam::symmTensor& A,
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
    const scalar tmp5 = A.xx()*A.xx();
    const scalar tmp6 = A.xy()*A.xy();
    const scalar tmp7 = A.xz()*A.xz();
    const scalar tmp8 = tmp6 + tmp7;
    const scalar tmp9 = tmp5 + tmp8;
    const scalar tmp10 = 0.5*A.xz();
    const scalar tmp11 = A.xx()*A.xz();
    const scalar tmp12 = A.xy()*A.yz();
    const scalar tmp13 = A.xz()*A.zz();
    const scalar tmp14 = tmp11 + tmp12 + tmp13;
    const scalar tmp15 = 0.5*A.xx();
    const scalar tmp16 = A.yy()*A.yy();
    const scalar tmp17 = A.yz()*A.yz();
    const scalar tmp18 = A.zz()*A.zz();
    const scalar tmp19 = A.xx()*A.yy();
    const scalar tmp20 = tmp19*A.zz();
    const scalar tmp21 = A.xz()*A.yz();
    const scalar tmp22 = tmp21*A.xy();
    const scalar tmp23 = pow(A.xx() + A.yy() + A.zz(), 2);
    const scalar tmp24 = 0.5*tmp16 + tmp17 + 0.5*tmp18 - 0.5*tmp23 + 0.5*tmp5 + tmp8;
    const scalar tmp25 = -tmp24;
    const scalar tmp26 = tmp25*tmp25;
    const scalar tmp27 = tmp25*tmp25*tmp25;
    const scalar tmp28 = tmp25*tmp25*tmp25*tmp25;
    const scalar tmp29 = tmp25*tmp25*tmp25*tmp25*tmp25;
    const scalar tmp30 = tmp17*A.xx();
    const scalar tmp31 = tmp6*A.zz();
    const scalar tmp32 = tmp7*A.yy();
    const scalar tmp33 = -tmp30 - tmp31 - tmp32 + A.xx()*A.yy()*A.zz() + 2.0*A.xy()*A.xz()*A.yz();
    const scalar tmp34 = tmp33*tmp33;
    const scalar tmp35 = tmp33*tmp33*tmp33;
    const scalar tmp36 = tmp33*tmp33*tmp33*tmp33;
    const scalar tmp37 = tmp33*tmp33*tmp33*tmp33*tmp33;
    const scalar tmp38 = tmp33*(66965.964947122498*tmp16 + 133931.929894245*tmp17 + 66965.964947122498*tmp18 - 66965.964947122498*tmp23 + 66965.964947122498*tmp5 + 133931.929894245*tmp6 + 133931.929894245*tmp7);
    const scalar tmp39 = 1068.7762439282301*tmp16 + 2137.5524878564602*tmp17 + 1068.7762439282301*tmp18 - 1068.7762439282301*tmp23 + 1068.7762439282301*tmp5 + 2137.5524878564602*tmp6 + 2137.5524878564602*tmp7;
    const scalar tmp40 = 1.4*tmp33*tmp39;
    const scalar tmp41 = tmp34*(4958153.453709905*tmp16 + 9916306.9074198101*tmp17 + 4958153.453709905*tmp18 - 4958153.453709905*tmp23 + 4958153.453709905*tmp5 + 9916306.9074198101*tmp6 + 9916306.9074198101*tmp7);
    const scalar tmp42 = 929746.52961154003*tmp16 + 1859493.0592230801*tmp17 + 929746.52961154003*tmp18 - 929746.52961154003*tmp23 + 929746.52961154003*tmp5 + 1859493.0592230801*tmp6 + 1859493.0592230801*tmp7;
    const scalar tmp43 = 1.4*tmp34*tmp42;
    const scalar tmp44 = tmp35*(16850541.013691049*tmp16 + 33701082.027382098*tmp17 + 16850541.013691049*tmp18 - 16850541.013691049*tmp23 + 16850541.013691049*tmp5 + 33701082.027382098*tmp6 + 33701082.027382098*tmp7);
    const scalar tmp45 = 7395951.3822101001*tmp16 + 14791902.7644202*tmp17 + 7395951.3822101001*tmp18 - 7395951.3822101001*tmp23 + 7395951.3822101001*tmp5 + 14791902.7644202*tmp6 + 14791902.7644202*tmp7;
    const scalar tmp46 = 1.4*tmp35*tmp45;
    const scalar tmp47 = -243546726.4462975*tmp16 - 487093452.89259499*tmp17 - 243546726.4462975*tmp18 + 243546726.4462975*tmp23 - 243546726.4462975*tmp5 - 487093452.89259499*tmp6 - 487093452.89259499*tmp7;
    const scalar tmp48 = 1.4*tmp36*tmp47;
    const scalar tmp49 = tmp36*(-785383511.86101997*tmp16 - 1570767023.7220399*tmp17 - 785383511.86101997*tmp18 + 785383511.86101997*tmp23 - 785383511.86101997*tmp5 - 1570767023.7220399*tmp6 - 1570767023.7220399*tmp7);
    const scalar tmp50 = tmp26*tmp33;
    const scalar tmp51 = tmp28*tmp33;
    const scalar tmp52 = tmp26*tmp34;
    const scalar tmp53 = tmp33*(30250.3056757796*tmp16 + 60500.611351559201*tmp17 + 30250.3056757796*tmp18 - 30250.3056757796*tmp23 + 30250.3056757796*tmp5 + 60500.611351559201*tmp6 + 60500.611351559201*tmp7);
    const scalar tmp54 = tmp34*(-2995332.4334491799*tmp16 - 5990664.8668983597*tmp17 - 2995332.4334491799*tmp18 + 2995332.4334491799*tmp23 - 2995332.4334491799*tmp5 - 5990664.8668983597*tmp6 - 5990664.8668983597*tmp7);
    const scalar tmp55 = tmp35*(-283712455.50391853*tmp16 - 567424911.00783706*tmp17 - 283712455.50391853*tmp18 + 283712455.50391853*tmp23 - 283712455.50391853*tmp5 - 567424911.00783706*tmp6 - 567424911.00783706*tmp7);
    const scalar tmp56 = tmp36*(-3962634249.4110899*tmp16 - 7925268498.8221798*tmp17 - 3962634249.4110899*tmp18 + 3962634249.4110899*tmp23 - 3962634249.4110899*tmp5 - 7925268498.8221798*tmp6 - 7925268498.8221798*tmp7);
    const scalar tmp57 = tmp27*tmp33;
    const scalar tmp58 = tmp27*tmp34;
    const scalar tmp59 = tmp26*tmp35;
    const scalar tmp60 = -157.59180543833327*tmp16 + 5843.1014121627168*tmp17*A.xx() - 315.18361087666653*tmp17 - 157.59180543833327*tmp18 - 5843.1014121627168*tmp20 - 11686.202824325434*tmp22 + 157.59180543833327*tmp23 + 6901136.4625648148*tmp26*tmp35 - 2431.652733354043*tmp26 + 1754664.0928393295*tmp27*tmp33 + 48583772.653129004*tmp27*tmp34 + 8607.3235851231639*tmp27 - 11204.663193945982*tmp28 - 760.13765306859932*tmp29 - 873078.28093423136*tmp34 - 9875809.5555594936*tmp35 + 294728034.84386063*tmp36 + 5408425687.0350599*tmp37 - 0.8*tmp38 - tmp40 - 0.8*tmp41 - tmp43 - 0.8*tmp44 - tmp46 - tmp48 - 0.8*tmp49 - 157.59180543833327*tmp5 - 698944.86700303352*tmp50 - 1368943.5740614603*tmp51 - 40573628.86741145*tmp52 + 5843.1014121627168*tmp6*A.zz() - 315.18361087666653*tmp6 + 5843.1014121627168*tmp7*A.yy() - 315.18361087666653*tmp7 - (0.66666666666666667*tmp16 + 1.3333333333333333*tmp17 + 0.66666666666666667*tmp18 - 0.66666666666666667*tmp23 + 0.66666666666666667*tmp5 + 1.3333333333333333*tmp6 + 1.3333333333333333*tmp7 + 1.0)*(247.22882602352038*tmp16 + 494.45765204704077*tmp17 + 247.22882602352038*tmp18 + 6879.1151319721794*tmp20 + 13758.230263944359*tmp22 - 247.22882602352038*tmp23 + 3834.6384047050678*tmp26 - 13278.832221188281*tmp27 + 15476.047042344479*tmp28 + 5601.210975515819*tmp29 - 6879.1151319721794*tmp30 - 6879.1151319721794*tmp31 - 6879.1151319721794*tmp32 - 58225.536367777553*tmp34 - 55265229.681683518*tmp35 - 1833034259.478168*tmp36 - 15366093393.53508*tmp37 + 247.22882602352038*tmp5 - 57260.848802108041*tmp50 - 5991860.9531104071*tmp51 - 66792738.808180198*tmp52 + 1.2*tmp53 + 1.2*tmp54 + 1.2*tmp55 + 1.2*tmp56 + 2436515.5238744877*tmp57 + 158549793.7719996*tmp58 - 1948319935.4517961*tmp59 + 494.45765204704077*tmp6 + 494.45765204704077*tmp7 + 28.097554988519878) - 19.256621626473745;
    const scalar tmp61 = 206.02402168626699*tmp16 + 412.04804337253398*tmp17 + 206.02402168626699*tmp18 + 5732.5959433101498*tmp20 + 11465.1918866203*tmp22 - 206.02402168626699*tmp23 + 3195.53200392089*tmp26 - 11065.6935176569*tmp27 + 12896.7058686204*tmp28 + 4667.6758129298496*tmp29 - 5732.5959433101498*tmp30 - 5732.5959433101498*tmp31 - 5732.5959433101498*tmp32 - 48521.280306481298*tmp34 - 46054358.0680696*tmp35 - 1527528549.56514*tmp36 - 12805077827.9459*tmp37 + 206.02402168626699*tmp5 - 47717.374001756703*tmp50 - 4993217.4609253397*tmp51 - 55660615.673483498*tmp52 + tmp53 + tmp54 + tmp55 + tmp56 + 2030429.6032287399*tmp57 + 132124828.143333*tmp58 - 1623599946.20983*tmp59 + 412.04804337253398*tmp6 + 412.04804337253398*tmp7 + 23.414629157099899;
    const scalar tmp62 = tmp61*tmp9;
    const scalar tmp63 = 217.55057658016449*tmp16 + 435.10115316032898*tmp17 + 217.55057658016449*tmp18 + 7034.4365791647597*tmp20 + 14068.873158329519*tmp22 - 217.55057658016449*tmp23 + 3723.89335663877*tmp26 - 15939.2396237307*tmp27 + 32221.941625641699*tmp28 - 23215.348852529802*tmp29 - 7034.4365791647597*tmp30 - 7034.4365791647597*tmp31 - 7034.4365791647597*tmp32 + 823995.18736610597*tmp34 + 8009700.2684979597*tmp35 - 257258805.87056699*tmp36 - 3957693983.0447302*tmp37 + tmp38 + tmp41 + tmp44 + tmp49 + 217.55057658016449*tmp5 + 880683.51532791601*tmp50 + 2144190.90344474*tmp51 + 37901059.935526699*tmp52 - 2370104.5868925201*tmp57 - 44927559.185148999*tmp58 - 21313392.0223355*tmp59 + 435.10115316032898*tmp6 + 435.10115316032898*tmp7 + 24.940908165785999;
    const scalar tmp64 = -11.749039875570251*tmp16 - 23.498079751140502*tmp17 - 11.749039875570251*tmp18 + 11.749039875570251*tmp23 + 7323414.9421357801*tmp26*tmp34 + 7249697.9680739902*tmp26*tmp35 - 391.04425139783802*tmp26 + 101013.983339062*tmp27*tmp33 + 2960.0486527581402*tmp27 - 10409.2072189767*tmp28 + 13808.8690964946*tmp29 - 153.96582059350601*tmp30 - 153.96582059350601*tmp31 - 153.96582059350601*tmp32 + tmp33*tmp39 + tmp34*tmp42 + 152772.950743819*tmp34 + tmp35*tmp45 + 2477178.1005436601*tmp35 + tmp36*tmp47 - 63514992.962433599*tmp36 - 1601621786.1423399*tmp37 - 11.749039875570251*tmp5 - 4001.3894709281199*tmp50 - 247435.10621023699*tmp51 - 9029803.7892927192*tmp58 - 23.498079751140502*tmp6 - 23.498079751140502*tmp7 + 153.96582059350601*A.xx()*A.yy()*A.zz() + 307.93164118701202*A.xy()*A.xz()*A.yz() - 0.49721779011075401;
    const scalar tmp65 = 43.510115316032902*tmp16 + 87.020230632065804*tmp17 + 43.510115316032902*tmp18 + 1406.887315832952*tmp20 + 2813.7746316659041*tmp22 - 43.510115316032902*tmp23 + 744.77867132775407*tmp26 - 3187.8479247461401*tmp27 + 6444.3883251283405*tmp28 - 4643.0697705059602*tmp29 - 1406.887315832952*tmp30 - 1406.887315832952*tmp31 - 1406.887315832952*tmp32 + 164799.03747322119*tmp34 + 1601940.0536995921*tmp35 - 51451761.1741134*tmp36 - 791538796.60894608*tmp37 + 0.2*tmp38 + 0.2*tmp41 + 0.2*tmp44 + 0.2*tmp49 + 43.510115316032902*tmp5 + 176136.70306558322*tmp50 + 428838.18068894802*tmp51 + 7580211.9871053398*tmp52 - 474020.91737850406*tmp57 - 8985511.8370297998*tmp58 - 4262678.4044671003*tmp59 + 87.020230632065804*tmp6 + 87.020230632065804*tmp7 + 4.9881816331572004;
    const scalar tmp66 = -0.85714285714285714*tmp61*(-0.4*tmp16 - 0.8*tmp17 - 0.4*tmp18 + 0.4*tmp23 - 1.6*tmp26 - 0.66666666666666667*tmp30 - 0.66666666666666667*tmp31 - 0.66666666666666667*tmp32 - 0.4*tmp5 - 0.8*tmp6 - 0.8*tmp7 + 0.66666666666666667*A.xx()*A.yy()*A.zz() + 1.3333333333333333*A.xy()*A.xz()*A.yz() - 0.2) - 0.85714285714285714*tmp65*(-2.0*tmp16 - 4.0*tmp17 - 2.0*tmp18 + 2.0*tmp23 - 2.0*tmp5 - 4.0*tmp6 - 4.0*tmp7 + 1.0) + 0.85714285714285714*(tmp24 + 0.16666666666666667)*(-16.448655825798351*tmp16 - 32.897311651596702*tmp17 - 16.448655825798351*tmp18 + 16.448655825798351*tmp23 + 10252780.918990092*tmp26*tmp34 + 10149577.155303586*tmp26*tmp35 - 547.46195195697317*tmp26 + 141419.5766746868*tmp27*tmp33 + 4144.0681138613963*tmp27 - 14572.89010656738*tmp28 + 19332.41673509244*tmp29 - 215.55214883090841*tmp30 - 215.55214883090841*tmp31 - 215.55214883090841*tmp32 + 213882.13104134658*tmp34 + 3468049.3407611242*tmp35 - 88920990.14740704*tmp36 - 2242270500.5992756*tmp37 + tmp40 + tmp43 + tmp46 + tmp48 - 16.448655825798351*tmp5 - 5601.9452592993675*tmp50 - 346409.14869433176*tmp51 - 12641725.305009807*tmp58 - 32.897311651596702*tmp6 - 32.897311651596702*tmp7 + 215.55214883090841*A.xx()*A.yy()*A.zz() + 431.10429766181682*A.xy()*A.xz()*A.yz() - 0.6961049061550556) + 0.85714285714285714;
    const scalar tmp67 = tmp10*tmp66 + tmp64*(0.5*tmp11 + 0.5*tmp12 + 0.5*tmp13);
    const scalar tmp68 = tmp11*tmp63 + tmp14*tmp62 + tmp60*(tmp10*tmp9 + tmp14*tmp15) + tmp67;
    const scalar tmp69 = eigenVector.xx()*eigenVector.xx();
    const scalar tmp70 = eigenVector.xy()*eigenVector.xy();
    const scalar tmp71 = eigenVector.yx()*eigenVector.yx();
    const scalar tmp72 = eigenVector.yy()*eigenVector.yy();
    const scalar tmp73 = eigenVector.zx()*eigenVector.zx();
    const scalar tmp74 = eigenVector.zy()*eigenVector.zy();
    const scalar tmp75 = tmp69*tmp70 + tmp71*tmp72 + tmp73*tmp74;
    const scalar tmp76 = 0.33333333333333333*A.xz();
    const scalar tmp77 = 0.16666666666666667*A.xx();
    const scalar tmp78 = 0.16666666666666667*A.zz();
    const scalar tmp79 = tmp16 + tmp17 + tmp6;
    const scalar tmp80 = 0.33333333333333333*tmp14;
    const scalar tmp81 = A.xy()*A.xz();
    const scalar tmp82 = A.yy()*A.yz();
    const scalar tmp83 = A.yz()*A.zz();
    const scalar tmp84 = tmp81 + tmp82 + tmp83;
    const scalar tmp85 = A.xx()*A.xy();
    const scalar tmp86 = A.xy()*A.yy();
    const scalar tmp87 = tmp21 + tmp85 + tmp86;
    const scalar tmp88 = 0.66666666666666667*tmp87;
    const scalar tmp89 = 0.16666666666666667*A.xz();
    const scalar tmp90 = 0.16666666666666667*A.yy();
    const scalar tmp91 = tmp60*(tmp14*tmp90 + tmp79*tmp89 + 0.33333333333333333*tmp84*A.xy() + 0.33333333333333333*tmp87*A.yz()) + tmp61*(tmp79*tmp80 + tmp84*tmp88) + tmp63*(0.66666666666666667*tmp12 + tmp76*A.yy()) + tmp64*(0.16666666666666667*tmp12 + tmp77*A.xz() + tmp78*A.xz()) + tmp66*tmp89;
    const scalar tmp92 = tmp75*tmp91;
    const scalar tmp93 = eigenVector.xz()*eigenVector.xz();
    const scalar tmp94 = eigenVector.yz()*eigenVector.yz();
    const scalar tmp95 = eigenVector.zz()*eigenVector.zz();
    const scalar tmp96 = tmp69*tmp93 + tmp71*tmp94 + tmp73*tmp95;
    const scalar tmp97 = tmp17 + tmp18 + tmp7;
    const scalar tmp98 = 0.5*A.zz();
    const scalar tmp99 = tmp61*tmp97;
    const scalar tmp100 = tmp13*tmp63 + tmp14*tmp99 + tmp60*(tmp10*tmp97 + tmp14*tmp98) + tmp67;
    const scalar tmp101 = tmp100*tmp96;
    const scalar tmp102 = eigenVector.xx()*eigenVector.xx()*eigenVector.xx();
    const scalar tmp103 = eigenVector.yx()*eigenVector.yx()*eigenVector.yx();
    const scalar tmp104 = eigenVector.zx()*eigenVector.zx()*eigenVector.zx();
    const scalar tmp105 = tmp102*eigenVector.xy() + tmp103*eigenVector.yy() + tmp104*eigenVector.zy();
    const scalar tmp106 = 0.33333333333333333*A.xx();
    const scalar tmp107 = 0.33333333333333333*tmp9;
    const scalar tmp108 = 0.16666666666666667*A.yz();
    const scalar tmp109 = tmp108*tmp66 + tmp60*(tmp108*tmp9 + tmp76*tmp87 + tmp77*tmp84 + tmp80*A.xy()) + tmp61*(tmp107*tmp84 + tmp14*tmp88) + tmp63*(tmp106*A.yz() + 0.66666666666666667*tmp81) + tmp64*(tmp78*A.yz() + 0.16666666666666667*tmp81 + tmp90*A.yz());
    const scalar tmp110 = tmp105*tmp109;
    const scalar tmp111 = tmp102*eigenVector.xz() + tmp103*eigenVector.yz() + tmp104*eigenVector.zz();
    const scalar tmp112 = 0.66666666666666667*tmp14;
    const scalar tmp113 = 0.16666666666666667*tmp17 + 0.16666666666666667*tmp5;
    const scalar tmp114 = 0.16666666666666667*tmp18 + 0.16666666666666667*tmp6;
    const scalar tmp115 = tmp64*(0.26666666666666667*tmp16 + 0.53333333333333333*tmp17 + 0.26666666666666667*tmp18 - 0.93333333333333333*tmp20 - 1.8666666666666667*tmp22 - 0.26666666666666667*tmp23 + 0.93333333333333333*tmp30 + 0.93333333333333333*tmp31 + 0.93333333333333333*tmp32 + 0.26666666666666667*tmp5 + 0.53333333333333333*tmp6 + 0.53333333333333333*tmp7 + 0.2);
    const scalar tmp116 = tmp61*(0.057142857142857143*tmp16 + 0.11428571428571429*tmp17 + 0.057142857142857143*tmp18 - 0.22857142857142857*tmp20 - 0.45714285714285714*tmp22 - 0.057142857142857143*tmp23 + 0.22857142857142857*tmp26 + 0.22857142857142857*tmp30 + 0.22857142857142857*tmp31 + 0.22857142857142857*tmp32 + tmp33*(-0.53333333333333333*tmp16 - 1.0666666666666667*tmp17 - 0.53333333333333333*tmp18 + 0.53333333333333333*tmp23 - 0.53333333333333333*tmp5 - 1.0666666666666667*tmp6 - 1.0666666666666667*tmp7) + 0.057142857142857143*tmp5 + 0.11428571428571429*tmp6 + 0.11428571428571429*tmp7 + 0.028571428571428571);
    const scalar tmp117 = tmp65*(-0.28571428571428571*tmp16 - 0.57142857142857143*tmp17 - 0.28571428571428571*tmp18 + 0.28571428571428571*tmp23 - 2.6666666666666667*tmp30 - 2.6666666666666667*tmp31 - 2.6666666666666667*tmp32 - 0.28571428571428571*tmp5 - 0.57142857142857143*tmp6 - 0.57142857142857143*tmp7 + 2.6666666666666667*A.xx()*A.yy()*A.zz() + 5.3333333333333333*A.xy()*A.xz()*A.yz() + 0.14285714285714286);
    const scalar tmp118 = -0.2*tmp115 - 0.2*tmp116 + 0.2*tmp117 - 0.028571428571428571;
    const scalar tmp119 = tmp118 + tmp60*(tmp112*A.xz() + tmp77*tmp97 + tmp78*tmp9) + tmp61*(tmp107*tmp97 + 0.66666666666666667*(tmp14*tmp14)) + tmp63*(tmp106*A.zz() + 0.66666666666666667*tmp7) + tmp64*(tmp113 + tmp114 + 0.33333333333333333*tmp7) + tmp66*(tmp77 + tmp78);
    const scalar tmp120 = tmp111*tmp119;
    const scalar tmp121 = 0.33333333333333333*A.zz();
    const scalar tmp122 = 0.33333333333333333*tmp97;
    const scalar tmp123 = 0.16666666666666667*A.xy();
    const scalar tmp124 = tmp123*tmp66 + tmp60*(tmp123*tmp97 + tmp76*tmp84 + tmp78*tmp87 + tmp80*A.yz()) + tmp61*(tmp112*tmp84 + tmp122*tmp87) + tmp63*(tmp121*A.xy() + 0.66666666666666667*tmp21) + tmp64*(0.16666666666666667*tmp21 + tmp77*A.xy() + tmp90*A.xy());
    const scalar tmp125 = tmp69*eigenVector.xy()*eigenVector.xz() + tmp71*eigenVector.yy()*eigenVector.yz() + tmp73*eigenVector.zy()*eigenVector.zz();
    const scalar tmp126 = 2.0*tmp125;
    const scalar tmp127 = tmp124*tmp126;
    const scalar tmp128 = -tmp102*eigenValue.x()*eigenVector.xz() - tmp103*eigenValue.y()*eigenVector.yz() - tmp104*eigenValue.z()*eigenVector.zz();
    const scalar tmp129 = 2.0*D.xz();
    const scalar tmp130 = 0.5*A.xy();
    const scalar tmp131 = tmp130*tmp66 + tmp64*(0.5*tmp21 + 0.5*tmp85 + 0.5*tmp86);
    const scalar tmp132 = tmp131 + tmp60*(tmp130*tmp9 + tmp15*tmp87) + tmp62*tmp87 + tmp63*tmp85;
    const scalar tmp133 = 0.5*A.yy();
    const scalar tmp134 = tmp61*tmp79;
    const scalar tmp135 = tmp131 + tmp134*tmp87 + tmp60*(tmp130*tmp79 + tmp133*tmp87) + tmp63*tmp86;
    const scalar tmp136 = tmp135*tmp75;
    const scalar tmp137 = tmp124*tmp96;
    const scalar tmp138 = 0.16666666666666667*tmp16 + 0.16666666666666667*tmp7;
    const scalar tmp139 = tmp118 + tmp60*(tmp77*tmp79 + tmp88*A.xy() + tmp9*tmp90) + tmp61*(tmp107*tmp79 + 0.66666666666666667*(tmp87*tmp87)) + tmp63*(0.33333333333333333*tmp19 + 0.66666666666666667*tmp6) + tmp64*(tmp113 + tmp138 + 0.33333333333333333*tmp6) + tmp66*(tmp77 + tmp90);
    const scalar tmp140 = tmp105*tmp139;
    const scalar tmp141 = tmp109*tmp111;
    const scalar tmp142 = tmp126*tmp91;
    const scalar tmp143 = -tmp102*eigenValue.x()*eigenVector.xy() - tmp103*eigenValue.y()*eigenVector.yy() - tmp104*eigenValue.z()*eigenVector.zy();
    const scalar tmp144 = 2.0*D.xy();
    const scalar tmp145 = 0.5*A.yz();
    const scalar tmp146 = tmp145*tmp66 + tmp64*(0.5*tmp81 + 0.5*tmp82 + 0.5*tmp83);
    const scalar tmp147 = tmp134*tmp84 + tmp146 + tmp60*(tmp133*tmp84 + tmp145*tmp79) + tmp63*tmp82;
    const scalar tmp148 = tmp146 + tmp60*(tmp145*tmp97 + tmp84*tmp98) + tmp63*tmp83 + tmp84*tmp99;
    const scalar tmp149 = 2.0*tmp105;
    const scalar tmp150 = 2.0*tmp111;
    const scalar tmp151 = tmp118 + tmp60*(tmp78*tmp79 + 0.66666666666666667*tmp84*A.yz() + tmp90*tmp97) + tmp61*(tmp122*tmp79 + 0.66666666666666667*(tmp84*tmp84)) + tmp63*(tmp121*A.yy() + 0.66666666666666667*tmp17) + tmp64*(tmp114 + tmp138 + 0.33333333333333333*tmp17) + tmp66*(tmp78 + tmp90);
    const scalar tmp152 = -tmp69*eigenValue.x()*eigenVector.xy()*eigenVector.xz() - tmp71*eigenValue.y()*eigenVector.yy()*eigenVector.yz() - tmp73*eigenValue.z()*eigenVector.zy()*eigenVector.zz();
    const scalar tmp153 = 2.0*D.yz();
    const scalar tmp154 = tmp105*tmp132;
    const scalar tmp155 = tmp111*tmp68;
    const scalar tmp156 = tmp109*tmp125;
    const scalar tmp157 = 2.0*tmp156;
    const scalar tmp158 = -0.6*tmp115 - 0.6*tmp116 + 0.6*tmp117 - 0.085714285714285714;
    const scalar tmp159 = tmp158 + tmp5*tmp63 + tmp60*tmp9*A.xx() + tmp61*(tmp9*tmp9) + tmp64*tmp9 + tmp66*A.xx();
    const scalar tmp160 = tmp139*tmp75;
    const scalar tmp161 = tmp119*tmp96;
    const scalar tmp162 = tmp158 + tmp18*tmp63 + tmp60*tmp97*A.zz() + tmp61*(tmp97*tmp97) + tmp64*tmp97 + tmp66*A.zz();
    const scalar tmp163 = -tmp69*tmp93*eigenValue.x() - tmp71*tmp94*eigenValue.y() - tmp73*tmp95*eigenValue.z();
    const scalar tmp164 = tmp158 + tmp16*tmp63 + tmp60*tmp79*A.yy() + tmp61*(tmp79*tmp79) + tmp64*tmp79 + tmp66*A.yy();
    const scalar tmp165 = -tmp69*tmp70*eigenValue.x() - tmp71*tmp72*eigenValue.y() - tmp73*tmp74*eigenValue.z();
    const scalar tmp166 = eigenVector.xy()*eigenVector.xy()*eigenVector.xy();
    const scalar tmp167 = eigenVector.yy()*eigenVector.yy()*eigenVector.yy();
    const scalar tmp168 = eigenVector.zy()*eigenVector.zy()*eigenVector.zy();
    const scalar tmp169 = tmp166*eigenVector.xx() + tmp167*eigenVector.yx() + tmp168*eigenVector.zx();
    const scalar tmp170 = tmp169*tmp91;
    const scalar tmp171 = tmp93*eigenVector.xx()*eigenVector.xy() + tmp94*eigenVector.yx()*eigenVector.yy() + tmp95*eigenVector.zx()*eigenVector.zy();
    const scalar tmp172 = tmp100*tmp171;
    const scalar tmp173 = tmp109*tmp75;
    const scalar tmp174 = tmp119*tmp125;
    const scalar tmp175 = tmp70*eigenVector.xx()*eigenVector.xz() + tmp72*eigenVector.yx()*eigenVector.yz() + tmp74*eigenVector.zx()*eigenVector.zz();
    const scalar tmp176 = 2.0*tmp175;
    const scalar tmp177 = tmp124*tmp176;
    const scalar tmp178 = tmp148*tmp171;
    const scalar tmp179 = tmp151*tmp175;
    const scalar tmp180 = -tmp70*eigenValue.x()*eigenVector.xx()*eigenVector.xz() - tmp72*eigenValue.y()*eigenVector.yx()*eigenVector.yz() - tmp74*eigenValue.z()*eigenVector.zx()*eigenVector.zz();
    const scalar tmp181 = tmp135*tmp169;
    const scalar tmp182 = tmp124*tmp171;
    const scalar tmp183 = tmp175*tmp91;
    const scalar tmp184 = 2.0*tmp183;
    const scalar tmp185 = tmp139*tmp169;
    const scalar tmp186 = tmp119*tmp171;
    const scalar tmp187 = tmp132*tmp75;
    const scalar tmp188 = tmp125*tmp68;
    const scalar tmp189 = tmp109*tmp176;
    const scalar tmp190 = tmp151*tmp171;
    const scalar tmp191 = tmp147*tmp175;
    const scalar tmp192 = -tmp166*eigenValue.x()*eigenVector.xx() - tmp167*eigenValue.y()*eigenVector.yx() - tmp168*eigenValue.z()*eigenVector.zx();
    const scalar tmp193 = 2.0*tmp124;
    const scalar tmp194 = -tmp93*eigenValue.x()*eigenVector.xx()*eigenVector.xy() - tmp94*eigenValue.y()*eigenVector.yx()*eigenVector.yy() - tmp95*eigenValue.z()*eigenVector.zx()*eigenVector.zy();
    const scalar tmp195 = tmp129*(tmp0*(-tmp105*tmp68 - tmp152 - tmp170 - tmp172 - 2.0*tmp173 - 2.0*tmp174 - tmp177) + tmp109) + tmp144*(tmp0*(-tmp154 - tmp157 - 2.0*tmp160 - tmp165 - tmp181 - tmp182 - tmp184) + tmp139) + tmp153*(tmp0*(-tmp110 - tmp127 - tmp147*tmp169 - tmp178 - 2.0*tmp179 - tmp180 - 2.0*tmp92) + tmp91) + (tmp0*(-2.0*tmp136 - tmp140 - tmp142 - tmp164*tmp169 - tmp190 - 2.0*tmp191 - tmp192) + tmp135)*D.yy() + (tmp0*(-tmp100*tmp126 - tmp105*tmp119 - tmp148*tmp176 - tmp151*tmp169 - tmp162*tmp171 - tmp193*tmp75 - tmp194) + tmp124)*D.zz() + (tmp0*(-tmp105*tmp159 - tmp143 - tmp185 - tmp186 - 2.0*tmp187 - 2.0*tmp188 - tmp189) + tmp132)*D.xx();
    const scalar tmp196 = eigenVector.xz()*eigenVector.xz()*eigenVector.xz();
    const scalar tmp197 = eigenVector.yz()*eigenVector.yz()*eigenVector.yz();
    const scalar tmp198 = eigenVector.zz()*eigenVector.zz()*eigenVector.zz();
    const scalar tmp199 = tmp196*eigenVector.xx() + tmp197*eigenVector.yx() + tmp198*eigenVector.zx();
    const scalar tmp200 = tmp124*tmp199;
    const scalar tmp201 = tmp135*tmp175;
    const scalar tmp202 = tmp109*tmp96;
    const scalar tmp203 = tmp125*tmp139;
    const scalar tmp204 = 2.0*tmp171;
    const scalar tmp205 = tmp204*tmp91;
    const scalar tmp206 = tmp100*tmp199;
    const scalar tmp207 = 2.0*tmp182;
    const scalar tmp208 = tmp119*tmp199;
    const scalar tmp209 = tmp139*tmp175;
    const scalar tmp210 = tmp68*tmp96;
    const scalar tmp211 = tmp125*tmp132;
    const scalar tmp212 = tmp109*tmp204;
    const scalar tmp213 = -tmp196*eigenValue.x()*eigenVector.xx() - tmp197*eigenValue.y()*eigenVector.yx() - tmp198*eigenValue.z()*eigenVector.zx();
    const scalar tmp214 = 2.0*tmp91;
    const scalar tmp215 = tmp129*(tmp0*(-tmp155 - tmp157 - 2.0*tmp161 - tmp163 - tmp183 - tmp206 - tmp207) + tmp119) + tmp144*(tmp0*(-tmp111*tmp132 - tmp152 - tmp200 - tmp201 - 2.0*tmp202 - 2.0*tmp203 - tmp205) + tmp109) + tmp153*(tmp0*(-2.0*tmp137 - tmp141 - tmp142 - tmp148*tmp199 - 2.0*tmp190 - tmp191 - tmp194) + tmp124) + (tmp0*(-2.0*tmp101 - tmp120 - tmp127 - tmp162*tmp199 - 2.0*tmp178 - tmp179 - tmp213) + tmp100)*D.zz() + (tmp0*(-tmp111*tmp139 - tmp126*tmp135 - tmp147*tmp204 - tmp151*tmp199 - tmp164*tmp175 - tmp180 - tmp214*tmp96) + tmp91)*D.yy() + (tmp0*(-tmp111*tmp159 - tmp128 - tmp208 - tmp209 - 2.0*tmp210 - 2.0*tmp211 - tmp212) + tmp68)*D.xx();
    const scalar tmp216 = eigenVector.xy()*eigenVector.xy()*eigenVector.xy()*eigenVector.xy();
    const scalar tmp217 = eigenVector.yy()*eigenVector.yy()*eigenVector.yy()*eigenVector.yy();
    const scalar tmp218 = eigenVector.zy()*eigenVector.zy()*eigenVector.zy()*eigenVector.zy();
    const scalar tmp219 = tmp216 + tmp217 + tmp218;
    const scalar tmp220 = tmp70*tmp93 + tmp72*tmp94 + tmp74*tmp95;
    const scalar tmp221 = tmp148*tmp220;
    const scalar tmp222 = tmp166*eigenVector.xz() + tmp167*eigenVector.yz() + tmp168*eigenVector.zz();
    const scalar tmp223 = tmp151*tmp222;
    const scalar tmp224 = -tmp166*eigenValue.x()*eigenVector.xz() - tmp167*eigenValue.y()*eigenVector.yz() - tmp168*eigenValue.z()*eigenVector.zz();
    const scalar tmp225 = tmp124*tmp220;
    const scalar tmp226 = tmp222*tmp91;
    const scalar tmp227 = 2.0*tmp169;
    const scalar tmp228 = tmp147*tmp222;
    const scalar tmp229 = tmp151*tmp220;
    const scalar tmp230 = 2.0*tmp222;
    const scalar tmp231 = -tmp70*tmp93*eigenValue.x() - tmp72*tmp94*eigenValue.y() - tmp74*tmp95*eigenValue.z();
    const scalar tmp232 = tmp196*eigenVector.xy() + tmp197*eigenVector.yy() + tmp198*eigenVector.zy();
    const scalar tmp233 = tmp124*tmp232;
    const scalar tmp234 = tmp220*tmp91;
    const scalar tmp235 = tmp148*tmp232;
    const scalar tmp236 = tmp151*tmp232;
    const scalar tmp237 = tmp147*tmp220;
    const scalar tmp238 = -tmp196*eigenValue.x()*eigenVector.xy() - tmp197*eigenValue.y()*eigenVector.yy() - tmp198*eigenValue.z()*eigenVector.zy();
    const scalar tmp239 = 2.0*tmp109;
    const scalar tmp240 = tmp129*(tmp0*(-tmp100*tmp232 - 2.0*tmp186 - tmp188 - tmp189 - tmp194 - 2.0*tmp225 - tmp226) + tmp124) + tmp144*(tmp0*(-tmp135*tmp222 - tmp180 - 2.0*tmp209 - tmp211 - tmp212 - tmp233 - 2.0*tmp234) + tmp91) + tmp153*(tmp0*(-tmp156 - tmp184 - tmp207 - tmp228 - 2.0*tmp229 - tmp231 - tmp235) + tmp151) + (tmp0*(-tmp119*tmp232 - tmp125*tmp159 - tmp132*tmp176 - tmp139*tmp222 - tmp152 - tmp204*tmp68 - tmp220*tmp239) + tmp109)*D.xx() + (tmp0*(-tmp162*tmp232 - 2.0*tmp172 - tmp174 - tmp177 - 2.0*tmp221 - tmp223 - tmp238) + tmp148)*D.zz() + (tmp0*(-tmp164*tmp222 - 2.0*tmp201 - tmp203 - tmp205 - tmp224 - tmp236 - 2.0*tmp237) + tmp147)*D.yy();
    const scalar tmp241 = eigenVector.xz()*eigenVector.xz()*eigenVector.xz()*eigenVector.xz();
    const scalar tmp242 = eigenVector.yz()*eigenVector.yz()*eigenVector.yz()*eigenVector.yz();
    const scalar tmp243 = eigenVector.zz()*eigenVector.zz()*eigenVector.zz()*eigenVector.zz();
    const scalar tmp244 = tmp241 + tmp242 + tmp243;


    result.xx() = tmp129*(tmp0*(-tmp101 - 2.0*tmp110 - 2.0*tmp120 - tmp127 - tmp128 - tmp4*tmp68 - tmp92) + tmp68) + tmp144*(tmp0*(-tmp132*tmp4 - tmp136 - tmp137 - 2.0*tmp140 - 2.0*tmp141 - tmp142 - tmp143) + tmp132) + tmp153*(tmp0*(-tmp109*tmp4 - tmp124*tmp150 - tmp126*tmp151 - tmp147*tmp75 - tmp148*tmp96 - tmp149*tmp91 - tmp152) + tmp109) + (tmp0*(-tmp100*tmp150 - tmp119*tmp4 - tmp124*tmp149 - tmp126*tmp148 - tmp151*tmp75 - tmp162*tmp96 - tmp163) + tmp119)*D.zz() + (tmp0*(-tmp126*tmp147 - tmp135*tmp149 - tmp139*tmp4 - tmp150*tmp91 - tmp151*tmp96 - tmp164*tmp75 - tmp165) + tmp139)*D.yy() + (tmp0*(tmp1*eigenValue.x() - 2.0*tmp154 - 2.0*tmp155 - tmp157 - tmp159*tmp4 - tmp160 - tmp161 + tmp2*eigenValue.y() + tmp3*eigenValue.z()) + tmp159)*D.xx();
    result.xy() = tmp195;
    result.xz() = tmp215;
    result.yy() = tmp129*(tmp0*(-tmp100*tmp220 - tmp109*tmp227 - tmp119*tmp176 - tmp180 - tmp193*tmp222 - tmp219*tmp91 - tmp68*tmp75) + tmp91) + tmp144*(tmp0*(-tmp135*tmp219 - 2.0*tmp185 - tmp187 - tmp189 - tmp192 - tmp225 - 2.0*tmp226) + tmp135) + tmp153*(tmp0*(-tmp147*tmp219 - 2.0*tmp170 - tmp173 - tmp177 - tmp221 - 2.0*tmp223 - tmp224) + tmp147) + (tmp0*(-tmp100*tmp176 - tmp119*tmp75 - tmp148*tmp230 - tmp151*tmp219 - tmp162*tmp220 - tmp169*tmp193 - tmp231) + tmp151)*D.zz() + (tmp0*(-tmp109*tmp230 - tmp119*tmp220 - tmp132*tmp227 - tmp139*tmp219 - tmp159*tmp75 - tmp165 - tmp176*tmp68) + tmp139)*D.xx() + (tmp0*(-tmp160 - tmp164*tmp219 - 2.0*tmp181 - tmp184 + tmp216*eigenValue.x() + tmp217*eigenValue.y() + tmp218*eigenValue.z() - 2.0*tmp228 - tmp229) + tmp164)*D.yy();
    result.yz() = tmp240;
    result.zz() = tmp129*(tmp0*(-tmp100*tmp244 - 2.0*tmp208 - tmp210 - tmp212 - tmp213 - 2.0*tmp233 - tmp234) + tmp100) + tmp144*(tmp0*(-tmp124*tmp244 - tmp132*tmp96 - tmp135*tmp220 - tmp139*tmp204 - tmp194 - tmp199*tmp239 - tmp214*tmp232) + tmp124) + tmp153*(tmp0*(-tmp148*tmp244 - 2.0*tmp200 - tmp202 - tmp205 - 2.0*tmp236 - tmp237 - tmp238) + tmp148) + (tmp0*(-tmp119*tmp244 - tmp132*tmp204 - tmp139*tmp220 - tmp159*tmp96 - tmp163 - 2.0*tmp199*tmp68 - tmp232*tmp239) + tmp119)*D.xx() + (tmp0*(-tmp135*tmp204 - tmp139*tmp96 - 2.0*tmp147*tmp232 - tmp151*tmp244 - tmp164*tmp220 - tmp199*tmp214 - tmp231) + tmp151)*D.yy() + (tmp0*(-tmp161 - tmp162*tmp244 - 2.0*tmp206 - tmp207 - tmp229 - 2.0*tmp235 + tmp241*eigenValue.x() + tmp242*eigenValue.y() + tmp243*eigenValue.z()) + tmp162)*D.zz();

}
