#ifndef CONFIGSETTINGS_H
#define CONFIGSETTINGS_H
#include <iostream>
#include <TMath.h> // needed for Long64_t type
using namespace std;

// ------------------------------------------------------
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
string cwd="/d/grid17/ln16/q-values-3";
string rootFileLoc="degALL_flatEtapi_b1_trees_subset_shap.root";
string rootTreeName="tree";
string fileTag="flatEtapi";
string runTag="";
string s_accWeight="AccWeight";
string s_sbWeight="weightBS";
string s_discrimVar="Mpi0;Meta";
string s_phaseVar="cosTheta_eta_gj;phi_eta_gj;cosTheta_X_cm;Phi;Mpi0eta;Mpi0g3";
string standardizationType="range";
string alwaysSaveTheseEvents="56815;63691;26045;21429;329;75269;44002;62359;9835;49405;83136;27544;60089;56015;33545;14328;87622;4280;45858;48721;90341;93993;76969;40429;32251;87630;93466;86693;88947;55028;81364;46104;3220;56577;91144;44551;86382;84076;20027;80983;38738;31799;22855;7927;79061;8367;7599;77494;33650;83600;622;76308;69750;43479;98437;98545;97394;55451;43255;96566;9362;18564;54645;446;71791;93101;41038;92326;4488;56925;46323;6916;55368;97191;42259;95145;96914;906;40039;88095;40263;81921;34190;65497;5591;92960;4923;36177;55423;76891;39016;3770;59274;66200;58721;85273;80420;94630;68902;27613;4417;74976;55522;19935;5180;5444;40827;65156;20654;97671;49175;30283;29137;99543;28748;44311;27558;24925;54953;23373;16354;18774;75361;30955;30071;50849;22508;97564;18517;18217;8006;7612;31043;40968;57638;46749;38823;40111;88831;19560;85497;66026;41347;65841;65093;59625;1635;94523;84351;60135;16392;32106;83423;94918;68375;53719;67349;76754;82163;32224;20226;78190;13086;19740;78249;45335;70433;57129;65692;17917;63080;65537;98423;629;13944;59916;68778;98081;64731;73716;93394;1017;22997;65031;52058;78250;22268;99833;62487;82391;51675;18712;5447;16861;19672;60293;67361;6224;67627;63373";
bool saveMemUsuage=1;
int nProcess=48;
int kDim=1000;
const int ckDim=1000; // same as kDim but just of const int type
const int phaseSpaceDim=6;
const int discrimVarDim=2;
bool redistributeBkgSigFits=0;
bool doKRandomNeighbors=0;
int numberEventsToSavePerProcess=0;
int seedShift=1341;
Long64_t nentries=-1;
int nRndRepSubset=0;
int nBS=0;
bool saveBShistsAlso=0;
bool override_nentries=0;
bool saveEventLevelProcessSpeed=1;
bool saveBranchOfNeighbors=1;
bool saveMemUsage=1;
bool verbose_outputDistCalc=false;
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
// ------------------------------------------------------
//
#endif
