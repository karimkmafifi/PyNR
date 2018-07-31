/*

Copyright (c) 2006-2010, The Scripps Research Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Author: Dr. Oleg Trott <ot14@columbia.edu>,
The Olson Lab,
The Scripps Research Institute

*/

#ifndef ATOM_CONSTANTS_H
#define ATOM_CONSTANTS_H

#include "bond.h"

//!Anzahl der internen Atomtypen:
const int n_intern_types = 160;

//!Versionsstring fuer die atom types:
const std::string a_t_version = "250911_2";

const int n_metal_list = 19;

const std::string metal_list[n_metal_list] = { // Si, Se and As not considered as metals
    "LI","BE","NA","MG","AL"," K","CA","CR","MN","FE","CO","NI","CU","ZN","AG","SN","CE","PT","AU"
};

const std::string metal_list_ra[n_metal_list] = { // Si, Se and As not considered as metals
    " LI"," BE"," NA"," MG"," AL","  K"," CA"," CR"," MN"," FE"," CO"," NI"," CU"," ZN"," AG"," SN"," CE"," PT"," AU"
};

//!es folgen die internen Atomtypen (mode 0):
const std::string i_t[n_intern_types] = {
    "H.ac","H.onh","H.n","H.o","H.0",                                //5

    "C.ar6p","C.ar6x","C.ar6","C.arp","C.arx","C.ar","C.2r3o","C.2r3x","C.3r3x","C.2r3","C.3r3",
    "C.1n","C.1p","C.1s","C.co2h","C.co2","C.es","C.hal","C.am","C.o","C.s","C.gu","C.guh",
    "C.mi","C.mih","C.n","C.2p","C.2s","C.2t","C.et","C.ohp","C.ohs","C.oht","C.3n","C.3p",
    "C.3s","C.3t","C.3q",                                        //38

    "N.ar6p","N.ar6","N.arp","N.ar2","N.ar3","N.ar3h","N.r3","N.az","N.1","N.o2","N.ohac","N.oh","N.ims",
    "N.imt","N.amp","N.ams","N.amt","N.samp","N.sams","N.samt","N.gu1","N.gu2","N.guh","N.mi1",
    "N.mi2","N.mih","N.aap","N.aas2","N.aas3","N.aat2","N.aat3","N.2n","N.2p","N.2s","N.2t","N.3n","N.3p","N.3s",
    "N.3t","N.4q","N.4h",                                            //41

    "O.ar","O.r3","O.h2o","O.n","O.noh","O.2co2","O.2es","O.2hal","O.am","O.co2","O.2po","O.2so",
    "O.2p","O.2s","O.3po","O.3so","O.carb","O.o","O.3ac","O.ph",
    "O.3oh","O.3es","O.3eta","O.3et",                                //24

    "S.ar","S.r3","S.thi","S.o","S.o2h","S.o3h","S.o4h","S.o2","S.o3",
    "S.o4","S.2","S.sh","S.s","S.3",                                //14

    "P.r3","P.o","P.o2h","P.o3h","P.o4h","P.o2","P.o3","P.o4","P.3",                //9

    "F.0","Cl.0","Br.0","I.0","F.i","Cl.i","Br.i","I.i",                        //8

    "Li","Na","Mg","Al","Si","K","Ca","Cr.th","Cr.oh","Mn","Fe","Co","Cu","Zn","Se","Mo","Sn",
    "Ni","Hg","B","As"                                            //21
};

//!es folgen die entsprechenden Sybyltypen (mode 1):
const std::string mode_1[n_intern_types] = {
    "H","H","H","H","H",                                        //5

    "C.ar","C.ar","C.ar","C.2","C.2","C.ar","C.2","C.2","C.3","C.2","C.3",
    "C.1","C.1","C.1","C.2","C.2","C.2","C.2","C.2","C.2","C.2","C.2","C.cat",
    "C.2","C.2","C.2","C.2","C.2","C.2","C.3","C.3","C.3","C.3","C.3","C.3",
    "C.3","C.3","C.3",                                        //38

    "N.pl3","N.ar","N.pl3","N.2","N.pl3","N.pl3","N.3","N.1","N.1","N.pl3","N.am","N.3","N.am",
    "N.am","N.am","N.am","N.am","N.am","N.am","N.am","N.2","N.pl3","N.pl3","N.2",
    "N.pl3","N.pl3","N.pl3","N.pl3","N.3","N.pl3","N.3","N.2","N.2","N.2","N.pl3","N.3","N.3","N.3",
    "N.3","N.4","N.4",                                            //41

    "O.3","O.3","O.3","O.2","O.3","O.2","O.2","O.2","O.2","O.co2","O.2","O.2",
    "O.co2","O.co2","O.3","O.3","O.2","O.3","O.3","O.3",
    "O.3","O.3","O.3","O.3",                                    //24

    "S.3","S.3","S.2","S.o","S.o2","S.o2","S.o2","S.o2","S.o2","S.o2","S.2","S.3","S.3","S.3",    //14

    "P.3","P.3","P.3","P.3","P.3","P.3","P.3","P.3","P.3",                        //9

    "F","Cl","Br","I","F","Cl","Br","I",                                //8

    "Li","Na","Mg","Al","Si","K","Ca","Cr.th","Cr.oh","Mn","Fe","Co","Cu","Zn","Se","Mo","Sn",
    "Ni","Hg","B","As"                                            //21
};

//!es folgen die modifizierten Sybyltypen (mode 2):
const std::string mode_2[n_intern_types] = {
    "H","H","H","H","H",                                        //5

    "C.ar","C.ar","C.ar","C.ar","C.ar","C.ar","C.2","C.2","C.3","C.2","C.3",
    "C.1","C.1","C.1","C.2","C.2","C.2","C.2","C.2","C.2","C.2","C.2","C.cat",
    "C.2","C.2","C.2","C.2","C.2","C.2","C.3","C.3","C.3","C.3","C.3","C.3",
    "C.3","C.3","C.3",                                        //38

    "N.pl3","N.ar","N.pl3","N.ar","N.pl3","N.pl3","N.3","N.1","N.1","N.pl3","N.am","N.3","N.am",
    "N.am","N.am","N.am","N.am","N.am","N.am","N.am","N.2","N.pl3","N.pl3","N.2",
    "N.pl3","N.pl3","N.pl3","N.pl3","N.3","N.pl3","N.3","N.2","N.2","N.2","N.pl3","N.3","N.3","N.3",
    "N.3","N.4","N.4",                                            //41

    "O.2","O.3","O.3","O.2","O.3","O.2","O.2","O.2","O.2","O.co2","O.2","O.2",
    "O.co2","O.co2","O.3","O.3","O.2","O.3","O.3","O.3",
    "O.3","O.3","O.3","O.3",                                    //24

    "S.2","S.3","S.2","S.o","S.o2","S.o2","S.o2","S.o2","S.o2","S.o2","S.2","S.3","S.3","S.3",    //14

    "P.3","P.3","P.3","P.3","P.3","P.3","P.3","P.3","P.3",                        //9

    "F","Cl","Br","I","F","Cl","Br","I",                                //8

    "Li","Na","Mg","Al","Si","K","Ca","Cr.th","Cr.oh","Mn","Fe","Co","Cu","Zn","Se","Mo","Sn",
    "Ni","Hg","B","As"                                            //21
};

//!Jetzt die Anzahl der Valenzen fuer jeden Atomtyp:
const int atom_valence[n_intern_types] = {
    1,1,1,1,1,                                            //5

    3,3,3,3,3,3,3,3,4,3,4,
    2,2,2,3,3,3,3,3,3,3,3,3,
    3,3,3,3,3,3,4,4,4,4,4,4,
    4,4,4,                                                //38

    3,2,3,2,3,3,3,2,1,3,3,3,3,
    3,3,3,3,3,3,3,2,3,3,2,
    3,3,3,3,3,3,3,2,2,2,3,3,3,3,
    3,4,4,                                                //41

    2,2,2,1,2,1,1,1,1,1,1,1,
    1,1,2,2,1,2,2,2,
    2,2,2,2,                                            //24

    2,2,1,3,4,4,4,4,4,
    4,1,2,2,2,                                            //14

    4,4,4,4,4,4,4,4,4,                                        //9

    1,1,1,1,0,0,0,0,                                        //8

    0,0,0,0,4,0,0,0,0,0,0,0,0,0,4,0,0,
    0,0,3,4                                                //21
};

//!Bindungsgeometrien:
const int atom_geometrie[n_intern_types] = { //!0:nichts / 1:linear / 2:trigonal planar / 3:tetraedrisch / 4:unknown
    0,0,0,0,0,                                            //5

    2,2,2,2,2,2,2,2,3,2,3,
    1,1,1,2,2,2,2,2,2,2,2,2,
    2,2,2,2,2,2,3,3,3,3,3,3,
    3,3,3,                                                //38

    2,2,2,2,2,2,3,1,1,2,2,3,2,
    2,2,2,2,2,2,2,2,3,2,2,
    3,2,3,2,3,2,3,2,2,2,2,3,3,3,
    3,3,3,                                                //41

    2,3,3,2,3,2,2,2,2,2,2,2,
    2,2,3,3,2,3,3,3,
    3,3,3,3,                                            //24

    2,3,3,2,3,3,3,3,3,
    3,2,3,3,3,                                            //14

    3,3,3,3,3,3,3,3,3,                                        //9

    0,0,0,0,0,0,0,0,                                        //8

    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4                                                //21
};

//!Hybridisierungen:
const int atom_hybridizations[n_intern_types] = {
    0,0,0,0,0,                                            //5

    2,2,2,2,2,2,2,2,3,2,3,
    1,1,1,2,2,2,2,2,2,2,2,2,
    2,2,2,2,2,2,3,3,3,3,3,3,
    3,3,3,                                                //38

    2,2,2,2,2,2,3,1,1,2,2,3,2,
    2,2,2,2,2,2,2,2,3,2,2,
    3,2,3,2,3,2,3,2,2,2,2,3,3,3,
    3,3,3,                                                //41

    2,3,3,2,3,2,2,2,2,2,2,2,
    2,2,3,3,2,3,3,3,
    3,3,3,3,                                            //24

    2,3,2,2,2,2,2,2,2,
    2,2,3,3,3,                                            //14

    3,2,2,2,2,2,2,2,3,                                        //9

    0,0,0,0,0,0,0,0,                                        //8

    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0                                                //21
};

//! X--H  Bindungslaengen:
const float atom_hb_length[n_intern_types] = { //!Typen die keine Protonen haben koennen haben dummy-werte!!!
    1.08061,1.08061,1.08061,1.08061,1.08061,                            //5

    0.956379,0.955857,0.960605,0.956379,0.952796,0.960605,0.968379,0.968379,0.979254,0.968379,0.976759,
    0.952907,0.952907,0.965793,0.991346,1.01374,0.955684,0.955684,0.982986,0.980125,0.978913,0.950231,0.950231,
    0.963721,0.950231,0.966764,0.972083,0.962652,0.962652,0.982888,0.989513,0.991182,0.991182,0.98106,0.975382,
    0.984575,0.985716,0.985716,                                    //38

    0.913131,0.913131,0.902918,0.899805,0.899805,0.899805,0.901294,0.897547,0.897547,0.909496,0.881489,0.909496,0.8869,
    0.8869,0.899901,0.886293,0.886293,0.865833,0.864477,0.864477,0.907462,0.887531,0.882072,0.884302,
    0.897229,0.889454,0.894858,0.889982,0.889982,0.889982,0.889982,0.899218,0.892104,0.902016,0.902016,0.902304,0.894173,0.889304,
    0.889304,0.922775,0.922775,                                        //41

    0.853137,1.10295,0.876543,0.914642,0.914642,0.92502,0.853137,0.853137,0.853137,0.92502,0.930701,0.949969,
    0.853137,0.752269,1.36144,1.36144,0.853137,0.917998,0.92502,0.892335,
    0.872473,0.92502,1.39311,1.10295,                                //24

    1.07966,1.23896,1.23896,1.23896,1.23896,1.23896,1.23896,1.23896,1.23896,
    1.23896,1.23896,1.2022,1.23896,1.23896,                                //14

    1.23896,1.23896,1.23896,1.23896,1.23896,1.23896,1.23896,1.23896,1.23896,            //9

    0.973445,1.24552,1.4,1.6,0.973445,1.24552,1.4,1.6,                        //8

    1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
    1.,1.,1.,1.                                            //21
};

const std::string normal_aacids[31] = { "ALA","ARG","ASN","ASP","ASH","ASX","CYS","CYX","CYM","GLN","GLU",
    "GLH","GLX","GLY","ILE","LEU","LYS","LYN","MET","PHE","PRO","SER",
    "THR","TRP","TYR","TYM","VAL","HIS","HID","HIE","HIP" };
const std::string modified_aacids[11] = { "ABA","CGU","CME","CSD","MEN","MLY","MSE","PCA","PTR","SEP","TPO" };

const std::string nucleic_acids[11] = { "  A","  T","  G","  C","  U"," +U"," DA"," DT"," DG"," DC"," DU" };

const std::string known_elements[29] = {
    "H","C","N","O","S","P","F","Cl","Br","I","Li","Na","Mg","Al","Si","K","Ca","Cr",
    "Mn","Fe","Co","Cu","Zn","Se","Mo","Sn","Ni","Hg","B"
};
const int number_of_elements = 25;
const std::string vdW_elements[number_of_elements] = {
    "H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I",
    "Li", "Na", "Mg", "Al", "Si", "K", "Ca", "B", "X", "Zn",
    "Fe", "Cu", "Ni",
    "Co", "Mn"
};
const float vdW_radii[number_of_elements] = { //united atom modell (wasserstoffe mit drin)
    1.20, //! H
    1.77, //! C
    1.35, //! N (um 0.2 erniedrigt um moegliche H-Bond zu beruecksichtigen)
    1.32, //! O (um 0.2 erniedrigt um moegliche H-Bond zu beruecksichtigen)
    1.65, //! S (um 0.2 erniedrigt um moegliche H-Bond zu beruecksichtigen)
    1.80, //! P
    1.47, //! F
    1.75, //! Cl
    1.85, //! Br
    1.98, //! I
    0.46, //! Li+   bei Koordinationszahl 6 (vdW - 0.3)
    0.72, //! Na+   bei Koordinationszahl 6 (vdW - 0.3)
    0.42, //! Mg++  bei Koordinationszahl 6 (vdW - 0.3)
    0.24, //! Al+++ bei Koordinationszahl 6 (vdW - 0.3)
    2.10, //! Si
    1.08, //! K+    bei Koordinationszahl 6 (vdW - 0.3)
    0.70, //! Ca++  bei Koordinationszahl 6 (vdW - 0.3)
    0.80, //! B
    1.60, //! X
    0.44, //! Zn++  bei Koordinationszahl 6 (vdW - 0.3)
    0.31, //! Fe++  bei Koordinationszahl 6 (vdW - 0.3)
    0.43, //! Cu++  bei Koordinationszahl 6 (vdW - 0.3)
    0.39, //! Ni++  bei Koordinationszahl 6 (vdW - 0.3)
    0.35, //! Co++  bei Koordinationszahl 6 (vdW - 0.3)
    0.53  //! Mn++  bei Koordinationszahl 6 (vdW - 0.3)
};

const float clash_radii[number_of_elements] = { //Atomradien, bzw. Ionenradien
                                                //! http://www.uniterra.de/rutherford/
    0.373, //! H
    0.772, //! C
    0.710, //! N
    0.604, //! O
    1.040, //! S
    0.930, //! P
    0.709, //! F
    0.994, //! Cl
    1.145, //! Br
    1.331, //! I
    0.780, //! Li+
    0.980, //! Na+
    0.780, //! Mg++
    0.570, //! Al+++
    1.170, //! Si
    1.330, //! K+
    1.060, //! Ca++
    0.830, //! B
    0.600, //! X
    0.830, //! Zn++
    0.670, //! Fe+++
    0.720, //! Cu++
    0.620, //! Ni+++
    0.640, //! Co++
    0.520  //! Mn+++
};
const float mweights[number_of_elements] = {
    1.008,  //! H
    12.011, //! C
    14.007, //! N
    15.999, //! O
    32.070, //! S
    30.974, //! P
    18.998, //! F
    35.453, //! Cl
    79.900, //! Br
    126.90, //! I
    6.941,  //! Li
    22.990, //! Na
    24.305, //! Mg
    26.982, //! Al
    28.086, //! Si
    39.100, //! K
    40.080, //! Ca
    10.810, //! B
    0.0,    //! X
    65.38,  //! Zn
    55.85,  //! Fe
    63.55,  //! Cu
    58.70,  //! Ni
    58.93,  //! Co
    54.93   //! Mn
};

// based on SY_TYPE_* but includes H
const int EL_TYPE_H = 0;
const int EL_TYPE_C = 1;
const int EL_TYPE_N = 2;
const int EL_TYPE_O = 3;
const int EL_TYPE_S = 4;
const int EL_TYPE_P = 5;
const int EL_TYPE_F = 6;
const int EL_TYPE_Cl = 7;
const int EL_TYPE_Br = 8;
const int EL_TYPE_I = 9;
const int EL_TYPE_Met = 10;
const int EL_TYPE_SIZE = 11;

// AutoDock4
const int AD_TYPE_C = 0;
const int AD_TYPE_A = 1;
const int AD_TYPE_N = 2;
const int AD_TYPE_O = 3;
const int AD_TYPE_P = 4;
const int AD_TYPE_S = 5;
const int AD_TYPE_H = 6; // non-polar hydrogen
const int AD_TYPE_F = 7;
const int AD_TYPE_I = 8;
const int AD_TYPE_NA = 9;
const int AD_TYPE_OA = 10;
const int AD_TYPE_SA = 11;
const int AD_TYPE_HD = 12;
const int AD_TYPE_Mg = 13;
const int AD_TYPE_Mn = 14;
const int AD_TYPE_Zn = 15;
const int AD_TYPE_Ca = 16;
const int AD_TYPE_Fe = 17;
const int AD_TYPE_Cl = 18;
const int AD_TYPE_Br = 19;
const int AD_TYPE_SIZE = 20;

const int AD_HBOND_TYPE_SIZE = 4;

// X-Score
const int XS_TYPE_C_H = 0;
const int XS_TYPE_C_P = 1;
const int XS_TYPE_N_P = 2;
const int XS_TYPE_N_D = 3;
const int XS_TYPE_N_A = 4;
const int XS_TYPE_N_DA = 5;
const int XS_TYPE_O_P = 6;
const int XS_TYPE_O_D = 7;
const int XS_TYPE_O_A = 8;
const int XS_TYPE_O_DA = 9;
const int XS_TYPE_S_A = 10;
const int XS_TYPE_S_P = 11;
const int XS_TYPE_P_P = 12;
const int XS_TYPE_F_H = 13;
const int XS_TYPE_Cl_H = 14;
const int XS_TYPE_Br_H = 15;
const int XS_TYPE_I_H = 16;
const int XS_TYPE_Met_D = 17;
const int XS_TYPE_SIZE = 18;


const int XS_HBOND_TYPE_SIZE = 7;

// Mol2 Atom Types
const int SY_TYPE_C_3 = 0;
const int SY_TYPE_C_2 = 1;
const int SY_TYPE_C_ar = 2;
const int SY_TYPE_C_1 = 3;
const int SY_TYPE_N_3 = 4;
const int SY_TYPE_N_2 = 5;
const int SY_TYPE_N_1 = 6;
const int SY_TYPE_O_3 = 7;
const int SY_TYPE_O_2 = 8;
const int SY_TYPE_S_3 = 9;
const int SY_TYPE_N_ar = 10;
const int SY_TYPE_P_3 = 11;
const int SY_TYPE_H = 12;
const int SY_TYPE_Br = 13;
const int SY_TYPE_Cl = 14;
const int SY_TYPE_F = 15;
const int SY_TYPE_I = 16;
const int SY_TYPE_S_2 = 17;
const int SY_TYPE_N_pl3 = 18;
const int SY_TYPE_LP = 19;
const int SY_TYPE_Na = 20;
const int SY_TYPE_K = 21;
const int SY_TYPE_Ca = 22;
const int SY_TYPE_Li = 23;
const int SY_TYPE_Al = 24;
const int SY_TYPE_Du = 25;
const int SY_TYPE_Du_C = 26;
const int SY_TYPE_Si = 27;
const int SY_TYPE_N_am = 28;
const int SY_TYPE_S_o = 29;
const int SY_TYPE_S_o2 = 30;
const int SY_TYPE_N_4 = 31;
const int SY_TYPE_O_co2 = 32;
const int SY_TYPE_C_cat = 33;
const int SY_TYPE_H_spc = 34;
const int SY_TYPE_O_spc = 35;
const int SY_TYPE_H_t3p = 36;
const int SY_TYPE_O_t3p = 37;
const int SY_TYPE_ANY = 38;
const int SY_TYPE_HEV = 39;
const int SY_TYPE_HET = 40;
const int SY_TYPE_HAL = 41;
const int SY_TYPE_Mg = 42;
const int SY_TYPE_Cr_oh = 43;
const int SY_TYPE_Cr_th = 44;
const int SY_TYPE_Se = 45;
const int SY_TYPE_Fe = 46;
const int SY_TYPE_Cu = 47;
const int SY_TYPE_Zn = 48;
const int SY_TYPE_Sn = 49;
const int SY_TYPE_Mo = 50;
const int SY_TYPE_Mn = 51;
const int SY_TYPE_Co_oh = 52;
const int SY_TYPE_SIZE = 53;


// Mol2 Bond Types
const int SY_BOND_1 = 0;
const int SY_BOND_2 = 1;
const int SY_BOND_3 = 2;
const int SY_BOND_am = 3;
const int SY_BOND_ar = 4;
const int SY_BOND_du = 5;
const int SY_BOND_un = 6;
const int SY_BOND_nc = 7;

struct Atom_kind {
    std::string name;
    float radius;
    float depth;
    float solvation;
    float volume;
    float covalent_radius; // from http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
    int hbond;
    float hbond_radius;
    float hbond_depth;
    int lone_pairs;
    float sasa_radii;
};

struct Metal_atom_kind_sasa {
    std::string name;
    float sasa_radii;
};

// generated from edited AD4_parameters.data using a script,
// then covalent radius added from en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
const Atom_kind atom_kind_data[] = { // name, radius, depth, solvation parameter, volume, covalent radius, hbond, hbond_radius, hbond_depth, lone_pairs, sasa_radii
    { "C",     2.00000,    0.15000,   -0.00143,   33.51030,    0.77,   0,   0,    0,    0,  1.70 }, //  0
    { "A",     2.00000,    0.15000,   -0.00052,   33.51030,    0.77,   0,   0,    0,    0,  1.70 }, //  1
    { "N",     1.75000,    0.16000,   -0.00162,   22.44930,    0.75,   0,   0,    0,    1,  1.55 }, //  2
    { "O",     1.60000,    0.20000,   -0.00251,   17.15730,    0.73,   0,   0,    0,    2,  1.52 }, //  3
    { "P",     2.10000,    0.20000,   -0.00110,   38.79240,    1.06,   0,   0,    0,    1,  1.80 }, //  4
    { "S",     2.00000,    0.20000,   -0.00214,   33.51030,    1.02,   0,   0,    0,    2,  1.80 }, //  5
    { "H",     1.00000,    0.02000,    0.00051,   0.00000,     0.37,   0,   0,    0,    0,  1.20 }, //  6
    { "F",     1.54500,    0.08000,   -0.00110,   15.44800,    0.71,   0,   0,    0,    3,  1.47 }, //  7
    { "I",     2.36000,    0.55000,   -0.00110,   55.05850,    1.33,   0,   0,    0,    3,  1.98 }, //  8
    { "NA",    1.75000,    0.16000,   -0.00162,   22.44930,    0.75,   4,   1.9,  5.0,  1,  1.55 }, //  9
    { "OA",    1.60000,    0.20000,   -0.00251,   17.15730,    0.73,   5,   1.9,  5.0,  2,  1.52 }, // 10
    { "SA",    2.00000,    0.20000,   -0.00214,   33.51030,    1.02,   5,   2.5,  1.0,  2,  1.80 }, // 11
    { "HD",    1.00000,    0.02000,    0.00051,   0.00000,     0.37,   2,   0,    0,    0,  1.20 }, // 12
    { "Mg",    0.65000,    0.87500,   -0.00110,   1.56000,     1.30,   0,   0,    0,    0,  1.73 }, // 13
    { "Mn",    0.65000,    0.87500,   -0.00110,   2.14000,     1.39,   0,   0,    0,    0,  1.73 }, // 14
    { "Zn",    0.74000,    0.55000,   -0.00110,   1.70000,     1.31,   0,   0,    0,    0,  1.39 }, // 15
    { "Ca",    0.99000,    0.55000,   -0.00110,   2.77000,     1.74,   0,   0,    0,    0,  1.80 }, // 16
    { "Fe",    0.65000,    0.01000,   -0.00110,   1.84000,     1.25,   0,   0,    0,    0,  1.80 }, // 17
    { "Cl",    2.04500,    0.27600,   -0.00110,   35.82350,    0.99,   0,   0,    0,    3,  1.75 }, // 18
    { "Br",    2.16500,    0.38900,   -0.00110,   42.56610,    1.14,   0,   0,    0,    3,  1.85 }  // 19
};

const Metal_atom_kind_sasa metal_atom_kind_sasa_data[] = { // name, sasa_radii
    { "Cu", 1.40 },
    { "Na", 2.27 },
    { "K",  2.75 },
    { "Hg", 1.80 }
};

const int metal_atom_kind_sasa_size = sizeof(metal_atom_kind_sasa_data) / sizeof(const Metal_atom_kind_sasa);

const float metal_solvation_parameter = -0.00110;

const float metal_covalent_radius = 1.75; // for metals not on the list // FIXME this info should be moved to non_ad_metals

const int atom_kinds_size = sizeof(atom_kind_data) / sizeof(const Atom_kind);

struct Atom_equivalence {
    std::string name;
    std::string to;
};

const Atom_equivalence atom_equivalence_data[] = {
    { "Se",  "S" }
};

const int atom_equivalences_size = sizeof(atom_equivalence_data) / sizeof(const Atom_equivalence);

const float xs_vdw_radii[] =
{
    1.9, // C_H
    1.9, // C_P
    1.8, // N_P
    1.8, // N_D
    1.8, // N_A
    1.8, // N_DA
    1.7, // O_P
    1.7, // O_D
    1.7, // O_A
    1.7, // O_DA
    2.0, // S_A
    2.0, // S_P
    2.1, // P_P
    1.5, // F_H
    1.8, // Cl_H
    2.0, // Br_H
    2.2, // I_H
    1.2  // Met_D
};

const std::string non_ad_metal_names[] = { // expand as necessary
    "Cu", "Na", "K", "Hg", "Co", "U", "Cd", "Ni"
};

const std::string xs_type_names[] = {
    "C_H", "C_P", "N_P", "N_D", "N_A", "N_DA", "O_P", "O_D", "O_A", "O_DA", "S_A", "S_P", "P_P", "F_H", "Cl_H", "Br_H", "I_H", "Met_D"
};

const int xs_kinds_size = sizeof(xs_type_names) / sizeof(const std::string);

static float max_covalent_rad = 0;

inline float metal_atom_sasa(const std::string& name) {

    for (int i = 0; i < metal_atom_kind_sasa_size; i++)
    {
        if (metal_atom_kind_sasa_data[i].name == name)
        {
            return metal_atom_kind_sasa_data[i].sasa_radii;
        }
    }

    return 1.80;
}

inline bool ad_is_hydrogen(int ad) {
    return ad == AD_TYPE_H || ad == AD_TYPE_HD;
}

inline bool ad_is_heteroatom(int ad) { // returns false for ad >= AD_TYPE_SIZE
    return ad != AD_TYPE_A && ad != AD_TYPE_C  &&
        ad != AD_TYPE_H && ad != AD_TYPE_HD &&
        ad < AD_TYPE_SIZE;
}

inline int ad_type_to_el_type(int t) {
    switch (t) {
    case AD_TYPE_C: return EL_TYPE_C;
    case AD_TYPE_A: return EL_TYPE_C;
    case AD_TYPE_N: return EL_TYPE_N;
    case AD_TYPE_O: return EL_TYPE_O;
    case AD_TYPE_P: return EL_TYPE_P;
    case AD_TYPE_S: return EL_TYPE_S;
    case AD_TYPE_H: return EL_TYPE_H;
    case AD_TYPE_F: return EL_TYPE_F;
    case AD_TYPE_I: return EL_TYPE_I;
    case AD_TYPE_NA: return EL_TYPE_N;
    case AD_TYPE_OA: return EL_TYPE_O;
    case AD_TYPE_SA: return EL_TYPE_S;
    case AD_TYPE_HD: return EL_TYPE_H;
    case AD_TYPE_Mg: return EL_TYPE_Met;
    case AD_TYPE_Mn: return EL_TYPE_Met;
    case AD_TYPE_Zn: return EL_TYPE_Met;
    case AD_TYPE_Ca: return EL_TYPE_Met;
    case AD_TYPE_Fe: return EL_TYPE_Met;
    case AD_TYPE_Cl: return EL_TYPE_Cl;
    case AD_TYPE_Br: return EL_TYPE_Br;
    case AD_TYPE_SIZE: return EL_TYPE_SIZE;
    default: Q_ASSERT(false);
    }
    return EL_TYPE_SIZE; // to placate the compiler in case of warnings - it should never get here though
}

inline float xs_radius(int t) {
    const int n = sizeof(xs_vdw_radii) / sizeof(const float);
    Q_ASSERT(n == XS_TYPE_SIZE);
    if (t == n)
    {
        return xs_vdw_radii[0];
    }
    Q_ASSERT(t < n);
    return xs_vdw_radii[t];
}

inline const std::string& xs_name(int i) {
    Q_ASSERT(XS_TYPE_SIZE == xs_kinds_size);
    Q_ASSERT(i < xs_kinds_size);
    return xs_type_names[i];
}

inline bool is_non_ad_metal_name(const std::string& name) {
    const int s = sizeof(non_ad_metal_names) / sizeof(const std::string);
    for (int i = 0; i < s; i++)
    {
        if (non_ad_metal_names[i] == name)
        {
            return true;
        }
    }
    return false;
}

inline bool xs_is_hydrophobic(int xs) {
    return xs == XS_TYPE_C_H ||
        xs == XS_TYPE_F_H ||
        xs == XS_TYPE_Cl_H ||
        xs == XS_TYPE_Br_H ||
        xs == XS_TYPE_I_H;
}

inline bool xs_is_Cl(int xs) {
    return xs == XS_TYPE_Cl_H;
}

inline bool xs_is_Br(int xs) {
    return xs == XS_TYPE_Br_H;
}

inline bool xs_is_I(int xs) {
    return xs == XS_TYPE_I_H;
}

inline bool xs_is_halogen(int xs) {
    return xs_is_Cl(xs) || xs_is_Br(xs) || xs_is_I(xs);
}

inline bool ad_is_acceptor(int ad) {
    return ad == AD_TYPE_NA ||
        ad == AD_TYPE_OA ||
        ad == AD_TYPE_SA;
}

inline bool ad_is_donor(int ad) {
    return ad == AD_TYPE_HD;
}

inline bool ad_donor_acceptor(int t1, int t2) {
    return ad_is_donor(t1) && ad_is_acceptor(t2);
}

inline bool ad_h_bond_possible(int t1, int t2) {
    return ad_donor_acceptor(t1, t2) || ad_donor_acceptor(t2, t1);
}

inline bool xs_is_acceptor(int xs) {
    return xs == XS_TYPE_N_A ||
        xs == XS_TYPE_N_DA ||
        xs == XS_TYPE_O_A ||
        xs == XS_TYPE_O_DA;
}

inline bool xs_is_donor(int xs) {
    return xs == XS_TYPE_N_D ||
        xs == XS_TYPE_N_DA ||
        xs == XS_TYPE_O_D ||
        xs == XS_TYPE_O_DA ||
        xs == XS_TYPE_Met_D;
}

inline bool xs_donor_acceptor(int t1, int t2) {
    return xs_is_donor(t1) && xs_is_acceptor(t2);
}

inline bool xs_h_bond_possible(int t1, int t2) {
    return xs_donor_acceptor(t1, t2) || xs_donor_acceptor(t2, t1);
}

inline bool xs_is_hal_acceptor(int xs) {
    return xs_is_acceptor(xs) || xs == XS_TYPE_S_A;
}

inline bool xs_hal_cl_bond_possible(int t1, int t2) {
    return (xs_is_Cl(t1) && xs_is_hal_acceptor(t2)) || (xs_is_Cl(t2) && xs_is_hal_acceptor(t1));
}

inline bool xs_hal_br_bond_possible(int t1, int t2) {
    return (xs_is_Br(t1) && xs_is_hal_acceptor(t2)) || (xs_is_Br(t2) && xs_is_hal_acceptor(t1));
}

inline bool xs_hal_i_bond_possible(int t1, int t2) {
    return (xs_is_I(t1) && xs_is_hal_acceptor(t2)) || (xs_is_I(t2) && xs_is_hal_acceptor(t1));
}

inline bool xs_hal_any_bond_possible(int t1, int t2) {
    return xs_hal_cl_bond_possible(t1, t2) || xs_hal_br_bond_possible(t1, t2) || xs_hal_i_bond_possible(t1, t2);
}

inline const Atom_kind& ad_type_property(int i) {
    Q_ASSERT(AD_TYPE_SIZE == atom_kinds_size);
    Q_ASSERT(i < atom_kinds_size);
    return atom_kind_data[i];
}

inline int string_to_ad_type(const std::string& name) { // returns AD_TYPE_SIZE if not found (no exceptions thrown, because metals unknown to AD4 are not exceptional)
    for (int i = 0; i < atom_kinds_size; i++)
    {
        if (atom_kind_data[i].name == name)
        {
            return i;
        }
    }

    for (int i = 0; i < atom_equivalences_size; i++)
    {
        if (atom_equivalence_data[i].name == name)
        {
            return string_to_ad_type(atom_equivalence_data[i].to);
        }
    }
    return AD_TYPE_SIZE;
}

inline float max_covalent_radius() {

    if (max_covalent_rad == 0)
    {
        for (int i = 0; i < atom_kinds_size; i++)
        {
            if (atom_kind_data[i].covalent_radius > max_covalent_rad)
            {
                max_covalent_rad = atom_kind_data[i].covalent_radius;
            }
        }
    }

    return max_covalent_rad;
}

#endif
