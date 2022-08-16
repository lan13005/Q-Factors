#ifndef CONFIGSETTINGS_H
#define CONFIGSETTINGS_H
#include <iostream>
#include <TMath.h> // needed for Long64_t type
using namespace std;

// ------------------------------------------------------
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
string cwd="/d/grid17/ln16/QFactors";
string rootFileLoc="flatetapi_b1_mc_081222/bkgndSample_recon_acc_flat_subset_1200k.root";
string rootTreeName="kin";
string fileTag="flatetapi_b1_11111";
string runTag="";
string s_fitWeight="AccWeight";
string s_sigWeight="AccWeight";
string s_altWeight="AccWeight;weightBS";
string s_discrimVar="Mpi0;Meta";
string s_extraVar="";
string s_neighborReqs="Meta>0.36;Meta<0.75;Mpi0>0.085;Mpi0<0.185";
string s_phaseVar="cosTheta_eta_gj;phi_eta_gj;Mpi0eta;Mpi0g3;Mpi0g4";
string s_circularVar="phi_eta_gj";
string standardizationType="range";
string alwaysSaveTheseEvents="";
bool saveMemUsuage=1;
int nProcess=48;
int kDim=200;
const int ckDim=200; // same as kDim but just of const int type
const int phaseSpaceDim=5;
const int discrimVarDim=2;
const int extraVarDim=0;
const int fitWeightsDim=1;
const int sigWeightsDim=1;
const int altWeightsDim=2;
bool redistributeBkgSigFits=0;
bool doKRandomNeighbors=0;
int numberEventsToSavePerProcess=2;
int seedShift=1341;
Long64_t nentries=-1;
int nRndRepSubset=0;
int nBS=0;
bool saveBShistsAlso=0;
bool override_nentries=0;
bool saveEventLevelProcessSpeed=1;
bool saveBranchOfNeighbors=0;
bool saveMemUsage=1;
bool verbose_outputDistCalc=false; // dont forget to set to false if running over a large dataset, there is a ton that is output if true
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
// ------------------------------------------------------
//
#endif
