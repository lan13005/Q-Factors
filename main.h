#ifndef MAIN_H
#define MAIN_H

#include "configSettings.h"
#include "configPDFs.h"
#include "auxilliary/helperFuncs.h"
#include "auxilliary/drawPlots/drawPlots.h"
#include <thread>
#include "TThread.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <TCanvas.h>
#include <TRandom.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TPolyLine3D.h>
#include <TLeaf.h>
#include "Math/MinimizerOptions.h"

TRandom rgen;
using namespace std;

class QFactorAnalysis{
	private:
		// variables to help load the data
	        TFile *dataFile;
		TTree *dataTree;
		Long64_t total_nentries;

                // Track memory usage
                ProcInfo_t pinfo;
                int startMemResident;
                int endMemResident;
	
		std::chrono::time_point<std::chrono::high_resolution_clock> start2;
                
                // These block of variables will be used to hold the initialization parameters. In the Q-factor paper they use 3 different initializations which
                // correspond to 100% bkg, 50/50, and 100% sig. If we want to do this here, the yields in the bkg and signal need to be modified. These vectors
                // will hold that information
                std::vector<float> sigFracs;

                // initialize vectors to hold the discriminating and phase space variables
                parseVarString parseDiscrimVars;
	        parseVarString parsePhaseSpace;
                parseVarString parseEventsToSave;
                std::vector<std::vector<float>> discrimVars; // holds all the data for discriminating variables
                std::vector<std::vector<float>> phaseSpaceVars; // holds all the data for the phase space varialbes

                float weight;
		std::vector<float> accWeights; 
                // Not all combinations will be a valid pairing. Suppose we only care about spectroscopically unique pairs, 
                // then we can fill phasePoint2PotentialNeighbor with only unique combos.
		std::vector<int> phasePoint2PotentialNeighbor; 
	
	public:
		QFactorAnalysis(){ start2 = std::chrono::high_resolution_clock::now(); };
		void initialize(string rootFileLoc, string rootTreeName);
		void loadFitParameters(string fitLocation,string cwd);
		void loadData();
		void runQFactorThreaded(int iProcess);

};

#endif
