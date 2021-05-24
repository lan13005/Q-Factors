#ifndef CONFIGSETTINGS_H
#define CONFIGSETTINGS_H
#include <iostream>
#include <TMath.h> // needed for Long64_t type
using namespace std;

// ------------------------------------------------------
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
string cwd="/d/grid17/ln16/q-values-3";
string rootFileLoc="degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat_pol000_090.root";
string rootTreeName="degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat";
string fileTag="000_090";
string runTag="_test";
string s_accWeight="AccWeight";
string s_sbWeight="weightBS";
string s_discrimVar="Mpi0;Meta";
string s_phaseVar="cosTheta_eta_gj;phi_eta_gj;cosTheta_X_cm;phi_X_lab";
string standardizationType="range";
string alwaysSaveTheseEvents="";
bool saveMemUsuage=1;
int nProcess=1;
int kDim=200;
const int ckDim=200; // same as kDim but just of const int type
const int phaseSpaceDim=4;
const int discrimVarDim=2;
bool redistributeBkgSigFits=0;
bool doKRandomNeighbors=0;
int numberEventsToSavePerProcess=2;
int seedShift=1341;
Long64_t nentries=300;
int nRndRepSubset=0;
int nBS=0;
bool saveBShistsAlso=0;
bool override_nentries=1;
bool saveEventLevelProcessSpeed=1;
bool saveBranchOfNeighbors=0;
bool saveMemUsage=1;
bool verbose_outputDistCalc=false;
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
// ------------------------------------------------------
//
#endif
