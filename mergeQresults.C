//The goal of this program is to read in all the q-value results and organizes them to be used when creating any histogram
#include "configSettings.h"
#include "auxilliary/helperFuncs.h"

void mergeQresults(){
        cout << "\n========================" << endl;
        cout << "Merging QFactors results with the original root file" << endl;
        cout << "1. Load qfactors result tree" << endl;
        cout << "2. Load input Data tree -- Clone -- Add new branches -- Fill" << endl;
        cout << "========================" << endl;

	// ---------------------------------------------------------------------------
	// ------------------------------Settting branch addresses for Q-Value Results
	// ---------------------------------------------------------------------------
	// Read in the qvalue data
	string qResultsFileName = "logs"+runTag+"/"+fileTag+"/resultsMerged_"+fileTag+".root";
	cout << "\n------------------\nOpening " << qResultsFileName << endl;
	TFile* qResultsFile = new TFile((qResultsFileName).c_str());
        TTree* qResultsTree;
	qResultsFile->GetObject(rootTreeName.c_str(),qResultsTree);
        cout << "Loaded the Q-factor results file with nentries = " << qResultsTree->GetEntries() << endl;

	float qvalue;
	int fitStatus;
	float conjugate_qvalue;
	float bestNLL;
	float worstNLL;
        float worst_qvalue;
        float qvalueBS_std;
        float eff_nentries;
        int kDim1; // kDim1 is an int whereas kDim is a const int so we can make an array out of it
        int neighbors[kDim];

        string s_discrimVar2 = replace_str(s_discrimVar,";","_");
	qResultsTree->SetBranchAddress(("qvalue_"+s_discrimVar2).c_str(),&qvalue);
	qResultsTree->SetBranchAddress(("fitStatus_"+s_discrimVar2).c_str(),&fitStatus);
	qResultsTree->SetBranchAddress(("bestNLL_"+s_discrimVar2).c_str(),&bestNLL);
	qResultsTree->SetBranchAddress(("worstNLL_"+s_discrimVar2).c_str(),&worstNLL);
        qResultsTree->SetBranchAddress(("worst_qvalue_"+s_discrimVar2).c_str(),&worst_qvalue);
        qResultsTree->SetBranchAddress(("qvalueBS_std_"+s_discrimVar2).c_str(),&qvalueBS_std);
        qResultsTree->SetBranchAddress(("eff_nentries_"+s_discrimVar2).c_str(),&eff_nentries);
        qResultsTree->SetBranchAddress(("kDim_"+s_discrimVar2).c_str(),&kDim1);
        if (saveBranchOfNeighbors)
            qResultsTree->SetBranchAddress(("neighbors_"+s_discrimVar2).c_str(),neighbors);
        cout << "Setting branch address..." << endl;

        int nentries;
	nentries=qResultsTree->GetEntries();
        cout << "Total number of events in Q-factor results: " << nentries << endl;

	// ---------------------------------------------------------------------------
	// ----------------------------------------------------- Loading the datafile data
	// and clone the datafile and add new branches to track the qvalue, NLL, etc
        // Now we have a root tree that has all the data from the DSelector and the q-factors analysis
	// ---------------------------------------------------------------------------
	string inputFileLoc = rootFileLoc;
	cout << "\n------------------\nOpening " << qResultsFileName << endl;
	TFile* dataFile=new TFile((inputFileLoc).c_str());
	TTree *dataTree;
	dataFile->GetObject((rootTreeName).c_str(),dataTree);
        cout << "Entries in input root tree: " << dataTree->GetEntries()  << endl;
	
	string postQFileName = "logs"+runTag+"/"+fileTag+"/postQVal_flatTree_"+fileTag+".root";	
	TFile *qd_dataFile = TFile::Open((postQFileName).c_str(),"RECREATE"); 
	TTree *outputTree = dataTree->CloneTree(-1,"fast"); 
	cout << "Cloned input to output file" << postQFileName << endl;

	TBranch* b_qvalue = outputTree->Branch(("qvalue_"+s_discrimVar2).c_str(),&qvalue,("qvalue_"+s_discrimVar2+"/F").c_str());
	TBranch* b_fitStatus = outputTree->Branch(("fitStatus_"+s_discrimVar2).c_str(),&fitStatus,("fitStatus_"+s_discrimVar2+"/I").c_str());
	TBranch* b_NLLBest = outputTree->Branch(("NLLBest_"+s_discrimVar2).c_str(),&bestNLL,("qvalue_NLLBest_"+s_discrimVar2+"/F").c_str());
	TBranch* b_NLLWorst = outputTree->Branch(("NLLWorst_"+s_discrimVar2).c_str(),&worstNLL,("qvalue_NLLWorst_"+s_discrimVar2+"/F").c_str());
	TBranch* b_worst_qvalue = outputTree->Branch(("worst_qvalue_"+s_discrimVar2).c_str(),&worst_qvalue,("worst_qvalue_"+s_discrimVar2+"/F").c_str());
	TBranch* b_qvalueBS_std = outputTree->Branch(("qvalueBS_std_"+s_discrimVar2).c_str(),&qvalueBS_std,("qvalueBS_std_"+s_discrimVar2+"/F").c_str());
        TBranch* b_eff_nentries = outputTree->Branch(("eff_nentries_"+s_discrimVar2).c_str(),&eff_nentries,("eff_nentries_"+s_discrimVar2+"/F").c_str());
        TBranch* b_kDim = outputTree->Branch(("kDim_"+s_discrimVar2).c_str(),&kDim1,("kDim_"+s_discrimVar2+"/I").c_str());
        TBranch* b_neighbors;
        if (saveBranchOfNeighbors)
            b_neighbors = outputTree->Branch(("neighbors_"+s_discrimVar2).c_str(),&neighbors,("neighbors_"+s_discrimVar2+"[kDim_"+s_discrimVar2+"]/I").c_str());
        cout << "Added new branches to output file..." <<endl;

	for (int ientry=0; ientry<nentries; ientry++)
	{
            qResultsTree->GetEntry(ientry);
	    b_qvalue->Fill();
	    b_fitStatus->Fill();
	    b_NLLBest->Fill();
	    b_NLLWorst->Fill();
	    b_worst_qvalue->Fill();
	    b_qvalueBS_std->Fill();
            b_eff_nentries->Fill();
            b_kDim->Fill();
            if (saveBranchOfNeighbors)
                b_neighbors->Fill();
	}
        cout << "Filled new output branches with qfactor results" << endl;
	qResultsFile->Close();
	qd_dataFile->cd();
	outputTree->Write(); 
        cout << "FINIHSED MERGING Q RESULTS!"<<endl;
}

