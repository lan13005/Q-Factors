#include "makePlots.h"

void makePlotsSummed(){
	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	
	TFile* dataFile2 = new TFile(("diagnosticPlots"+runTag+"/postQVal_hists.root").c_str());

        // Match histogram type to a given index
	TH1F* hists1D[5];
	TH2F* hists2D[5];
        map<string,int> types={
            {"_tot",0},
            {"_sig",1},
            {"_bkg",2},
            {"_sb_sig",3},
            {"_sb_bkg",4}
        };

        // FIRST LOOP THROUGH THE KEYS TO DETERMINE THE HISTOGRAMS TO PLOT
        set<string> names;
        std::string delimiter = "_";
        TIter nextkey(gDirectory->GetListOfKeys());
        while (TKey* key = (TKey*)nextkey() ){
            string className=key->ReadObj()->ClassName();
            string objName=key->GetName();
            cout << "Reading key: " << objName << " of class: " << className << endl;
            int tagLength=7;
            string tag = objName.substr(objName.size()-tagLength,objName.size());
            if (types.find(tag) == types.end()){
                tagLength=4;
                tag = objName.substr(objName.size()-tagLength,objName.size());
            }
            string name = objName.substr(0, objName.size()-tagLength);
            names.insert(name);
        }

        bool is1D=true;
        for (string nameToPlot : names){
            cout << nameToPlot << endl;
            TIter nextkey(gDirectory->GetListOfKeys());
            while (TKey* key = (TKey*)nextkey() ){
                string className=key->ReadObj()->ClassName();
                string objName=key->GetName();

                int tagLength=7;
                string tag = objName.substr(objName.size()-tagLength,objName.size());
                if (types.find(tag) == types.end()){
                    tagLength=4;
                    tag = objName.substr(objName.size()-tagLength,objName.size());
                }
                string name = objName.substr(0, objName.size()-tagLength);
                string type = objName.substr(objName.size()-tagLength,objName.length());

                bool condition=name.compare(nameToPlot)==0;
                if (condition){
                    cout << "type: " << type << endl;
                    if (className=="TH1F"){
                        cout << "  found match ("+objName+")! Histogram is TH1F" << endl;
                        is1D=true;
                        for (auto aType: types){
                            if (type==aType.first)
                                hists1D[types[type]]=(TH1F*)key->ReadObj();
                        }
                    }
                    else if (className=="TH2F"){
                        is1D=false;
                        cout << "  found match ("+objName+")! Histogram is TH2F" << endl;
                        for (auto aType: types){
                            if (type==aType.first)
                                hists2D[types[type]]=(TH2F*)key->ReadObj();
                        }
                    }
                    else { cout << "Unexpected key! Found object that is not TH1F nor TH2F. exiting..." << endl; exit(0); }
                }
            }
            if(is1D){
                cout << "Drawing hist: " << hists1D[0]->GetName() << endl;
	        makeStackedHist(hists1D[0],hists1D[1],hists1D[2],hists1D[3],hists1D[4],nameToPlot,"diagnosticPlots"+runTag);
            }
            else{
                cout << "Drawing hist: " << hists2D[0]->GetName() << endl;
                make2DHistsOnPads(hists2D[0],hists2D[1],hists2D[2],hists2D[3],hists2D[4],nameToPlot, "diagnosticPlots"+runTag);
            }
        }
}

