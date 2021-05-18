#include "makePlots.h"

void makePlotsSummed(){
	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	
	TFile* dataFile2 = new TFile(("diagnosticPlots"+runTag+"/postQVal_hists.root").c_str());

        // FIRST LOOP THROUGH THE KEYS TO DETERMINE THE HISTOGRAMS TO PLOT
        set<string> names;
        std::string delimiter = "_";
        TIter nextkey(gDirectory->GetListOfKeys());
        while (TKey* key = (TKey*)nextkey() ){
            string className=key->ReadObj()->ClassName();
            string objName=key->GetName();
            cout << "Reading key: " << objName << " of class: " << className << endl;
            string name = objName.substr(0, objName.find(delimiter));
            string type = objName.substr(objName.find(delimiter),objName.length());
            names.insert(name);
        }

	TH1F* hists1D[5];
	TH2F* hists2D[5];
        bool is1D=true;
        vector<string> types={"_tot","_sig","_bkg","_sig_sb","_bkg_sb"};
        for (string nameToPlot : names){
            TIter nextkey(gDirectory->GetListOfKeys());
            while (TKey* key = (TKey*)nextkey() ){
                string className=key->ReadObj()->ClassName();
                string objName=key->GetName();
                string name = objName.substr(0, objName.find(delimiter));
                if (name==nameToPlot){
                    string type = objName.substr(objName.find(delimiter),objName.length());
                    if (className=="TH1F"){
                        is1D=true;
                        for (int i=0; i<(int)types.size();++i){ 
                            if (type==types[i])
                                hists1D[i]=(TH1F*)key->ReadObj();
                        }
                    }
                    if (className=="TH2F"){
                        is1D=false;
                        for (int i=0; i<(int)types.size();++i){ 
                            if (type==types[i])
                                hists2D[i]=(TH2F*)key->ReadObj();
                        }
                    }
                }
            }
            if(is1D)
	        makeStackedHist(hists1D[0],hists1D[1],hists1D[2],hists1D[3],hists1D[4],nameToPlot,"diagnosticPlots"+runTag);
            else
                make2DHistsOnPads(hists2D[0],hists2D[1],hists2D[2],hists2D[3],hists2D[4],nameToPlot, "diagnosticPlots"+runTag);
        }
}

