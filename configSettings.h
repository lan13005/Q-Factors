#ifndef CONFIGSETTINGS_H
#define CONFIGSETTINGS_H
#include <iostream>
#include <TMath.h> // needed for Long64_t type
using namespace std;

// ------------------------------------------------------
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
string cwd="/d/grid17/ln16/q-values-3";
string rootFileLoc="/d/grid17/ln16/q-values-3/logs/flatEtapi_b1_Mpi0_1111/postQVal_flatTree_flatEtapi_b1_Mpi0_1111.root";
string rootTreeName="tree_4g_flat";
string fileTag="flatEtapi_b1_Mpi0_MetaIndep_1111";
string runTag="";
string s_fitWeight="AccWeight";
string s_sigWeight="AccWeight;qvalue_Mpi0";
string s_altWeight="AccWeight;weightBS";
string s_discrimVar="Meta";
string s_extraVar="";
string s_neighborReqs="";
string s_phaseVar="cosTheta_eta_gj;phi_eta_gj;Mpi0eta;Mpi0g3";
string standardizationType="range";
string alwaysSaveTheseEvents="";
bool saveMemUsuage=1;
int nProcess=2;
int kDim=400;
const int ckDim=400; // same as kDim but just of const int type
const int phaseSpaceDim=4;
const int discrimVarDim=1;
const int extraVarDim=0;
const int fitWeightsDim=1;
const int sigWeightsDim=2;
const int altWeightsDim=2;
bool redistributeBkgSigFits=0;
bool doKRandomNeighbors=0;
int numberEventsToSavePerProcess=2;
int seedShift=1341;
Long64_t nentries=500;
int nRndRepSubset=0;
int nBS=0;
bool saveBShistsAlso=0;
bool override_nentries=1;
bool saveEventLevelProcessSpeed=1;
bool saveBranchOfNeighbors=1;
bool saveMemUsage=1;
bool verbose_outputDistCalc=false;
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
// ------------------------------------------------------
//
#endif
