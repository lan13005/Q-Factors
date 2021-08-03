// This program loads data from the postQ root file and draws a bunch of histograms
#include "makePlots.h"
#include "configSettings.h"
#include "./auxilliary/helperFuncs.h"

//////////////////////
//////   INPUT  //////
//////////////////////
//Enter the branches to be plotted into this vector. For 2D histograms, separate the x and y variable by a semicolon
vector<string> histsToMake={
    "Meta",
    "Mpi0",
    "Mpi0g3",
    "Mpi0g4",
    "Mpi0eta",
    "cosTheta_X_cm",
    "cosTheta_eta_gj",
    "phi_eta_gj",
    "Mpi0;Meta",
    "Mpi0eta;cosTheta_eta_gj"
};

void makePlots(bool makeTotal){
	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
        gStyle->SetOptStat(0);

        ///////////////////////////////////////
        // LOAD THE DATA FROM THE POST-Q ROOT FILE
        ///////////////////////////////////////
	string postQFileName = "logs"+runTag+"/"+fileTag+"/postQVal_flatTree_"+fileTag+".root";	
        if (makeTotal)
            postQFileName = "logs"+runTag+"/"+"/postQVal_flatTree.root";
        cout << "Looking in " << postQFileName << endl;
	TFile* postQFile = TFile::Open((postQFileName).c_str()); 
        TTree* dataTree;
	postQFile->GetObject(rootTreeName.c_str(),dataTree);

        // SETUP VARIALBES TO LOAD QFACTOR RELATED VARIABLES
        cout << "LOADING Q-FATOR RELATED BRANCHES" << endl;
	if (!override_nentries){
	    nentries=dataTree->GetEntries();
            cout << "nentries: " << nentries << endl;
	}
	float qvalue;
	float bestNLL;
	float worstNLL;
        float worst_qvalue;
        float qvalueBS_std;
        float eff_nentries;
        bool is_truecombo;
        int neighbors[ckDim];
	dataTree->SetBranchAddress("qvalue",&qvalue);
	dataTree->SetBranchAddress("NLLBest",&bestNLL);
	dataTree->SetBranchAddress("NLLWorst",&worstNLL);
	dataTree->SetBranchAddress("worst_qvalue",&worst_qvalue);
	dataTree->SetBranchAddress("qvalueBS_std",&qvalueBS_std);
        dataTree->SetBranchAddress("eff_nentries",&eff_nentries);
        dataTree->SetBranchAddress("neighbors",&neighbors);
	std::vector< float > qvalues;
	std::vector< float > bestNLLs;
	std::vector< float > worstNLLs;
	std::vector< float > worst_qvalues;
	std::vector< float > qvalueBS_stds;
	std::vector< float > eff_nentrieses;
        std::vector< bool > is_truecombos;
        std::vector< std::array<int,ckDim> > neighborses;
        std::array<int,ckDim> copyableNeighbors;

        // SETUP VARIALBES FOR WEIGHTING
        cout << "LOADING WEIGHTING BRANCHES" << endl;
        double sbWeight;
        double accWeight;
        float sbWeight_f;
        float accWeight_f;
        Long64_t sbWeight_l;
        Long64_t accWeight_l;
        std::vector<float> accWeights; accWeights.reserve(nentries); 
	std::vector<float> sbWeights; sbWeights.reserve(nentries);
        string sbWeightType;
        string accWeightType;

        if(!s_accWeight.empty()){
            cout << "Using accidental weight branch" << endl;
            accWeightType=setBranchAddress(dataTree, s_accWeight, &accWeight_l, &accWeight_f, &accWeight);
            //dataTree->SetBranchAddress(s_accWeight.c_str(),&accWeight);
            //typeName=dataTree->GetLeaf(s_accWeight.c_str())->GetTypeName();
            //if (typeName=="Float_t")
            //    dataTree->SetBranchAddress(s_accWeight.c_str(),&accWeight_f);
            //if (typeName=="Double_t")
            //    dataTree->SetBranchAddress(s_accWeight.c_str(),&accWeight);
        }
        else{
            cout << "Not using accidental weights" << endl;
            accWeight_f=1;
            accWeightType="Float_t";
        }
        if(!s_sbWeight.empty()){
            cout << "Using sideband weight branch" << endl;
            sbWeightType=setBranchAddress(dataTree, s_sbWeight, &sbWeight_l, &sbWeight_f, &sbWeight);
            //dataTree->SetBranchAddress(s_sbWeight.c_str(),&sbWeight);
        }
        else{
            cout << "Not using sideband weights" << endl;
            sbWeight=1;
        }
        TBranch* br = (TBranch*)dataTree->GetListOfBranches()->FindObject("is_truecombo");
        if (br){
            dataTree->SetBranchAddress("is_truecombo",&is_truecombo);
            cout << "is_truecombo branch exists - will overlay the matchedThrown distributions for you" << endl;
        }
        
        // we dont want to load the branch name multiple times (i.e. if there are repeated names in histsToMake)
        // so we can insert them into a set while parsing the inputs to get information on what variables to plot and dimensionality (i.e. TH1 or TH2)
        cout << "DETERMINING WHICH VARIABLES TO PLOT AND DIMENSIONALITY OF HISTOGRAMS" << endl;
        vector<string> branchesToGet={};
        vector<vector<string>> varsToPlot;
        for (auto s: histsToMake){
            size_t pos = 0;
            string delim=";";
            std::string token;
            vector<string> tmp;
            while ((pos = s.find(delim)) != std::string::npos) {
                token = s.substr(0, pos);
                if (std::find(branchesToGet.begin(), branchesToGet.end(), token) == branchesToGet.end()) 
                    branchesToGet.push_back(token);
                tmp.push_back(token);
                s.erase(0, pos + delim.length());
            }
            tmp.push_back(s);
            varsToPlot.push_back(tmp);
            if (std::find(branchesToGet.begin(), branchesToGet.end(), s) == branchesToGet.end()) 
                branchesToGet.push_back(s);
        }
        cout << "LOADING THE BRANCHES TO PLOT: "<< endl;
        for ( auto ele : branchesToGet)
            cout << ele << endl;
        string typeName;
        vector<string> typeNames;
        vector<double> value(branchesToGet.size(),0);
        vector<float> value_f(branchesToGet.size(),0);
        vector<Long64_t> value_l(branchesToGet.size(),0);
        vector<float> minValue(branchesToGet.size(),FLT_MAX);
        vector<float> maxValue(branchesToGet.size(),FLT_MIN);
        vector<vector<float>> values;
        map<string,int> nameToIdx;
        int i=0;
        for (auto s: branchesToGet){
            //dataTree->SetBranchAddress(s.c_str(),&value[i]);
            typeName=setBranchAddress(dataTree, s, &value_l[i], &value_f[i], &value[i]);
            typeNames.push_back(typeName);
            nameToIdx[s]=i;
            values.push_back(vector<float>{});
            values[i].reserve(nentries);
            ++i;
        }

        cout << "LOADING DATA INTO MEMORY" << endl;
	for (int ientry=0; ientry<nentries; ientry++)
	{
		dataTree->GetEntry(ientry);

                if (accWeightType=="Float_t")
                    accWeights.push_back(accWeight_f);
                if (accWeightType=="Double_t")
                    accWeights.push_back(accWeight);
                if (accWeightType=="Long64_t")
                    accWeights.push_back(accWeight_l);
                if (sbWeightType=="Float_t")
		    sbWeights.push_back(sbWeight_f);
                if (sbWeightType=="Double_t")
		    sbWeights.push_back(sbWeight);
                if (sbWeightType=="Long64_t")
		    sbWeights.push_back(sbWeight_l);

                auto it = std::find(branchesToGet.begin(), branchesToGet.end(), "qvalue");
                if ( it  != branchesToGet.end() ){
                    //cout << "qvalue branch found in branchesToGet (set of variables used to make the requested histograms)" << endl;
                    int index=distance(branchesToGet.begin(),it); // index of the q-factors array
		    qvalues.push_back(value_f[index]); // qvalues is defined to be a float from the q-factors program
                }
                else{
                    qvalues.push_back(qvalue);
                }

                if (br)
                    is_truecombos.push_back(is_truecombo);


		bestNLLs.push_back(bestNLL);
		worstNLLs.push_back(worstNLL);
                worst_qvalues.push_back(worst_qvalue);
                qvalueBS_stds.push_back(qvalueBS_std);
                eff_nentrieses.push_back(eff_nentries);
                std::copy(std::begin(copyableNeighbors),std::end(copyableNeighbors),std::begin(neighbors));
                neighborses.push_back(copyableNeighbors);

                for (int j=0; j<(int)branchesToGet.size(); ++j){
                    if(typeNames[j]=="Float_t"){
                        values[j].push_back(value_f[j]);
                    }
                    if(typeNames[j]=="Double_t"){
                        value_f[j]=value[j]; // copy it over to value_f so we can simply code when finding min and max values below 
                        values[j].push_back(value_f[j]);
                    }
                    if(typeNames[j]=="Long64_t"){
                        value_f[j]=value_l[j]; 
                        values[j].push_back(value_f[j]);
                    }
                    if (value_f[j] < minValue[j])
                        minValue[j]=value_f[j];
                    if (value_f[j] > maxValue[j])
                        maxValue[j]=value_f[j];
                }
	}


        // -----------------------------------
        // Defining Histograms
        // -----------------------------------
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);

        cout << "INITIALIZING HISTOGRAMS" << endl;
        vector<TH1F*> hists1D_truth;
        vector<TH1F*> hists1D_tot;
        vector<TH1F*> hists1D_sig;
        vector<TH1F*> hists1D_bkg;
        vector<TH1F*> hists1D_sig_sb;
        vector<TH1F*> hists1D_bkg_sb;
        vector<TH2F*> hists2D_truth;
        vector<TH2F*> hists2D_tot;
        vector<TH2F*> hists2D_sig;
        vector<TH2F*> hists2D_bkg;
        vector<TH2F*> hists2D_sig_sb;
        vector<TH2F*> hists2D_bkg_sb;

        map<string, int> mapVarsToPlotToIdx;
        int idx1D=0;
        int idx2D=0;
        for (auto varToPlot: varsToPlot){
            if (varToPlot.size()==1){
                float minVal = minValue[nameToIdx[varToPlot[0]]];
                float maxVal = maxValue[nameToIdx[varToPlot[0]]];
                minVal=0.1;
                maxVal=3.0;
                hists1D_truth.push_back(new TH1F((varToPlot[0]+"_truth").c_str(),("truth;"+varToPlot[0]).c_str(),75,minVal,maxVal));
                hists1D_tot.push_back(new TH1F((varToPlot[0]+"_tot").c_str(),("tot;"+varToPlot[0]).c_str(),75,minVal,maxVal));
                hists1D_sig.push_back(new TH1F((varToPlot[0]+"_sig").c_str(),("sig;"+varToPlot[0]).c_str(),75,minVal,maxVal));
                hists1D_bkg.push_back(new TH1F((varToPlot[0]+"_bkg").c_str(),("bkg;"+varToPlot[0]).c_str(),75,minVal,maxVal));
                hists1D_sig_sb.push_back(new TH1F((varToPlot[0]+"_sb_sig").c_str(),("sb_sig;"+varToPlot[0]).c_str(),75,minVal,maxVal));
                hists1D_bkg_sb.push_back(new TH1F((varToPlot[0]+"_sb_bkg").c_str(),("sb_bkg;"+varToPlot[0]).c_str(),75,minVal,maxVal));
                mapVarsToPlotToIdx[varToPlot[0]]=idx1D;
                cout << "Setting " << varToPlot[0] << " to match with " << idx1D << endl;
                ++idx1D;
            }
            else if (varToPlot.size()==2){
                float minVal1 = minValue[nameToIdx[varToPlot[0]]];
                float maxVal1 = maxValue[nameToIdx[varToPlot[0]]];
                float minVal2 = minValue[nameToIdx[varToPlot[1]]];
                float maxVal2 = maxValue[nameToIdx[varToPlot[1]]];
                hists2D_truth.push_back(
                        new TH2F((varToPlot[0]+"Vs"+varToPlot[1]+"_truth").c_str(),("truth;"+varToPlot[0]+";"+varToPlot[1]).c_str(),75,minVal1,maxVal1,75,minVal2,maxVal2));
                hists2D_tot.push_back(
                        new TH2F((varToPlot[0]+"Vs"+varToPlot[1]+"_tot").c_str(),("tot;"+varToPlot[0]+";"+varToPlot[1]).c_str(),75,minVal1,maxVal1,75,minVal2,maxVal2));
                hists2D_sig.push_back(
                        new TH2F((varToPlot[0]+"Vs"+varToPlot[1]+"_sig").c_str(),("sig;"+varToPlot[0]+";"+varToPlot[1]).c_str(),75,minVal1,maxVal1,75,minVal2,maxVal2));
                hists2D_bkg.push_back(
                        new TH2F((varToPlot[0]+"Vs"+varToPlot[1]+"_bkg").c_str(),("bkg;"+varToPlot[0]+";"+varToPlot[1]).c_str(),75,minVal1,maxVal1,75,minVal2,maxVal2));
                hists2D_sig_sb.push_back(
                        new TH2F((varToPlot[0]+"Vs"+varToPlot[1]+"_sb_sig").c_str(),("sb_sig;"+varToPlot[0]+";"+varToPlot[1]).c_str(),
                            75,minVal1,maxVal1,75,minVal2,maxVal2));
                hists2D_bkg_sb.push_back(
                        new TH2F((varToPlot[0]+"Vs"+varToPlot[1]+"_sb_bkg").c_str(),("sb_bkg;"+varToPlot[0]+";"+varToPlot[1]).c_str(),
                            75,minVal1,maxVal1,75,minVal2,maxVal2));
                mapVarsToPlotToIdx[varToPlot[0]+";"+varToPlot[1]]=idx2D;
                cout << "Setting " << varToPlot[0]+";"+varToPlot[1] << " to match with " << idx2D << endl;
                ++idx2D;
            }
        }

        // We can now fill the histograms, properly weighted
        cout << "FILLING HISTOGRAMS" << endl;
	float sigWeight;
	float totWeight;
	float bkgWeight;
	float sigWeight_sb;
	float bkgWeight_sb;
        float baseWeight; 

	for (int ientry=0; ientry<nentries; ientry++){
		qvalue = qvalues[ientry];
                accWeight=accWeights[ientry];
                sbWeight=sbWeights[ientry];
                if (br)
                    is_truecombo=is_truecombos[ientry];
                baseWeight=accWeight;

                if (isnan(qvalue)){
                    cout << "(ientry:"<<ientry<<") qvalue is nan... exiting..." << endl;
                    exit(0);
                }

                ////////////////////////////////////
                // Multiply q-factor and accidetal weights if requested
                ////////////////////////////////////
		sigWeight = qvalue*baseWeight;
		totWeight = baseWeight;
		bkgWeight = totWeight-sigWeight;

                //////////////////////////////
                // Multiply sideband and accidental weights if requested 
                //////////////////////////////
                sigWeight_sb = baseWeight*sbWeight;
                bkgWeight_sb = totWeight-sigWeight_sb;
            
                ////////////////////////////////////
                // Fill histograms
                ////////////////////////////////////
                for (auto varToPlot: varsToPlot){
                    if (varToPlot.size()==1){
                        i=mapVarsToPlotToIdx[varToPlot[0]];
                        float val1=values[nameToIdx[varToPlot[0]]][ientry];
                        if ( (br!=NULL) & is_truecombo) 
                            hists1D_truth[i]->Fill(val1);
                        hists1D_tot[i]->Fill(val1,totWeight);
                        hists1D_sig[i]->Fill(val1,sigWeight);
                        hists1D_bkg[i]->Fill(val1,bkgWeight);
                        hists1D_sig_sb[i]->Fill(val1,sigWeight_sb);
                        hists1D_bkg_sb[i]->Fill(val1,bkgWeight_sb);
                    }
                    else if (varToPlot.size()==2){
                        i=mapVarsToPlotToIdx[varToPlot[0]+";"+varToPlot[1]];
                        float val1=values[nameToIdx[varToPlot[0]]][ientry];
                        float val2=values[nameToIdx[varToPlot[1]]][ientry];
                        if ( (br!=NULL) & is_truecombo) 
                            hists2D_truth[i]->Fill(val1,val2);
                        hists2D_tot[i]->Fill(val1,val2,totWeight);
                        hists2D_sig[i]->Fill(val1,val2,sigWeight);
                        hists2D_bkg[i]->Fill(val1,val2,bkgWeight);
                        hists2D_sig_sb[i]->Fill(val1,val2,sigWeight_sb);
                        hists2D_bkg_sb[i]->Fill(val1,val2,bkgWeight_sb);
                    }
                }
        }

	// HERE WE WILL JUST DRAW SOME OF THE HISTOGRAMS WITH THE BKG FILLED IN TO SEE THEIR CONTRIBUTION
        cout << "STACKING HISTOGRAMS INTO CANVAS" << endl;
        for (auto varToPlot: varsToPlot){
            if (varToPlot.size()==1){
                cout << "plotting " << varToPlot[0] << " with hist name: " << hists1D_tot[i]->GetName() << endl;
                i=mapVarsToPlotToIdx[varToPlot[0]];
                if (makeTotal)
	            makeStackedHist(hists1D_truth[i], hists1D_tot[i],hists1D_sig[i],hists1D_bkg[i],
                            hists1D_sig_sb[i],hists1D_bkg_sb[i],varToPlot[0], "diagnosticPlots"+runTag+"/");
                else 
	            makeStackedHist(hists1D_truth[i], hists1D_tot[i],hists1D_sig[i],hists1D_bkg[i],
                            hists1D_sig_sb[i],hists1D_bkg_sb[i],varToPlot[0], "diagnosticPlots"+runTag+"/"+fileTag);
            }
            else if (varToPlot.size()==2){
                i=mapVarsToPlotToIdx[varToPlot[0]+";"+varToPlot[1]];
                if (makeTotal)
                    make2DHistsOnPads(hists2D_truth[i], hists2D_tot[i],hists2D_sig[i],hists2D_bkg[i],hists2D_sig_sb[i],hists2D_bkg_sb[i],
                            varToPlot[0]+"Vs"+varToPlot[1], "diagnosticPlots"+runTag);
                else
                    make2DHistsOnPads(hists2D_truth[i], hists2D_tot[i],hists2D_sig[i],hists2D_bkg[i],hists2D_sig_sb[i],hists2D_bkg_sb[i],
                            varToPlot[0]+"Vs"+varToPlot[1], "diagnosticPlots"+runTag+"/"+fileTag);
            }
        }

        // we can also directly save all the histograms to a root file
        cout << "SAVING ALL HISTOGRAMS TO A ROOT FILE"<<endl;
        if (makeTotal)
	    TFile* dataFile3 = new TFile(("diagnosticPlots"+runTag+"/postQVal_hists.root").c_str(),"RECREATE");
        else
	    TFile* dataFile3 = new TFile(("diagnosticPlots"+runTag+"/"+fileTag+"/postQVal_hists_"+fileTag+".root").c_str(),"RECREATE");
        for (auto hist: hists1D_tot)
            hist->Write();
        for (auto hist: hists1D_sig)
            hist->Write();
        for (auto hist: hists1D_bkg)
            hist->Write();
        for (auto hist: hists1D_sig_sb)
            hist->Write();
        for (auto hist: hists1D_bkg_sb)
            hist->Write();
        for (auto hist: hists2D_tot)
            hist->Write();
        for (auto hist: hists2D_sig)
            hist->Write();
        for (auto hist: hists2D_bkg)
            hist->Write();
        for (auto hist: hists2D_sig_sb)
            hist->Write();
        for (auto hist: hists2D_bkg_sb)
            hist->Write();
}

