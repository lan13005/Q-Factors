#ifndef GETSBWEIGHT_H
#define GETSBWEIGHT_H

#include <iostream>
#include "./auxilliary/drawPlots/drawPlots.h"
#include "configSettings.h"
#include "TCanvas.h"
#include <RooAbsReal.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGenericPdf.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooConstVar.h>
#include <RooArgList.h>
#include <RooPlot.h>
#include <RooMinuit.h>
#include <RooFitResult.h>
#include <RooChi2Var.h>
#include <RooWorkspace.h>
#include <RooMinimizer.h>
#include "RooTrace.h"
//#include "RooAddPdf.h"
//#include "RooFormulaVar.h"

using namespace std;

class fitManager
{
    //////////////////////////////////////////////////////////////////////
    // 1. Need to define the values to initialize the fits with and the RooFit variables/PDFs
    // 2. Need to define a way to reinitialize the PDFs
    // 3. Need to define a way to calculate the q-factor
    // 4. Need to define how you want to insert the data into the dataset
    // 5. (DEFAULTS PROBABLY FINE) Need to define how you want to fit the PDFs to the dataset
    //////////////////////////////////////////////////////////////////////
    public:
        // DEFINE VALUES TO INITIALIZE ROOFIT VARIABLES WITH
        float initMassY = 0.549354;
        float initSigmaY = 0.0166494;
        float initBernC = 0.5;
        float initBernD = 0.5;
        std::vector<float> fitRangeY={0.36,0.75};

        ///////////////////////////////
        // DEFINE VARIABLES FOR DATASET
        ///////////////////////////////
        float eff_nentries;
        RooRealVar* roo_Meta;
        RooRealVar* roo_Weight;
        RooDataSet* rooData;
        ////////////////////////
        // DEFINE THE SIGNAL PDF
        ////////////////////////
        RooRealVar* peak_eta;
        RooRealVar* width_eta;
        RooGaussian* rooSig;
        ////////////////////////
        // DEFINE THE BKG PDF
        ////////////////////////
        RooRealVar* bern_parC;
        RooRealVar* bern_parD;
        RooGenericPdf* rooBkg;
        ////////////////////////
        // DEFINE THE TOTAL PDF
        ////////////////////////
        RooRealVar* nsig;
        RooRealVar* nbkg;
        RooAddPdf* rooSigPlusBkg;
        
	fitManager(string iProcess){ 
            /////////////// VARIABLES
            roo_Meta = new RooRealVar{("roo_Meta_"+iProcess).c_str(),"Mass GeV",(fitRangeY[1]-fitRangeY[0])/2+fitRangeY[0],fitRangeY[0],fitRangeY[1]};
            roo_Meta->setRange(("roo_fitRangeMeta_"+iProcess).c_str(),fitRangeY[0], fitRangeY[1]);
            roo_Meta->setBins(100);
            roo_Weight = new RooRealVar{("roo_Weight_"+iProcess).c_str(), "Weight", 0, -10, 10}; // Weights can take a wide range
            rooData = new RooDataSet{("rooData"+iProcess).c_str(),"rooData",RooArgSet(*roo_Meta,*roo_Weight),RooFit::WeightVar(*roo_Weight)};
            /////////////// FOR SIGNAL PDF
            peak_eta = new RooRealVar{("peak_eta_"+iProcess).c_str(),"peak_eta",initMassY};
            width_eta = new RooRealVar{("width_eta_"+iProcess).c_str(),"width_eta",initSigmaY,initSigmaY*0.5,initSigmaY*1.5};
            /////////////// FOR BKG PDF
            bern_parC = new RooRealVar{("bern_parC_"+iProcess).c_str(),"bern_parC",0.5,0,1};
            bern_parD = new RooRealVar{("bern_parD_"+iProcess).c_str(),"bern_parD",0.5,0,1};
            /////////////// FOR YIELDS
            nsig = new RooRealVar{("nsig_"+iProcess).c_str(),"nsig",(float)kDim/2,0,(float)kDim};
            nbkg = new RooRealVar{("nbkg_"+iProcess).c_str(),"nbkg",(float)kDim/2,0,(float)kDim};
            /////////////// CREATING FINAL PDFS
            // THE NAMES rooSig, rooBkg, rooSigPlusBkg CANNOT CHANGE
            rooSig = new RooGaussian{("rooSig_"+iProcess).c_str(), "rooGaus", *roo_Meta, *peak_eta, *width_eta};
            rooBkg = new RooGenericPdf{("rooBkg_"+iProcess).c_str(), "rooBkg"
                ,("bern_parC_"+iProcess+"*roo_Meta_"+iProcess+"+bern_parD_"+iProcess+"*(1-roo_Meta_"+iProcess+")").c_str()
                ,RooArgSet(*bern_parC,*bern_parD,*roo_Meta)};
            rooSigPlusBkg = new RooAddPdf{("rooSumPdf_"+iProcess).c_str(), "rooSumPdf", RooArgList(*rooSig,*rooBkg),RooArgSet(*nsig,*nbkg)};
            // THE NAMES rooSig, rooBkg, rooSigPlusBkg CANNOT CHANGE
            //////////////
        }

        void reinitialize(float initSigFrac){
            peak_eta->setVal(initMassY);
            width_eta->setVal(initSigmaY);
            bern_parC->setVal(initBernC);
            bern_parD->setVal(initBernD);
            eff_nentries = rooData->sumEntries();
            nsig->setVal(initSigFrac*eff_nentries);
            nbkg->setVal((1-initSigFrac)*eff_nentries);
        }

        float calculate_q(vector<vector<float>> &vect, int neighborIdx){
            roo_Meta->setVal(vect[0][neighborIdx]);
            float sigFrac = nsig->getVal()/(nsig->getVal()+nbkg->getVal());
            float sigPdfVal = sigFrac*rooSig->getVal(RooArgSet(*roo_Meta));
            float bkgPdfVal = (1-sigFrac)*rooBkg->getVal(RooArgSet(*roo_Meta));
            float totPdfVal = rooSigPlusBkg->getVal(RooArgSet(*roo_Meta));
            float qvalue = sigPdfVal/(sigPdfVal+bkgPdfVal);
            cout << "postFit(Q=" << qvalue << ") - nsig: " << nsig->getVal() << " || nbkg: " << nbkg->getVal() << endl;
            return qvalue;
        }

        bool insert(vector<vector<float>> &vect, int neighborIdx, float weight){
            // with 600 neighbors it seems like the fitTo command takes ~2x longer when using Range() argument which selects the fit range.
            // This is equivalent to shrinking the dataset range and fitting over the full range which will save time.
            // It might be useful to think of setting a fit range as to lower the effective number of nearest neighbors. Another thing that lowers
            // the effective number of neighbors is any weights we apply to the filling of the histograms
            float y = vect[0][neighborIdx];
            bool keptNeighbor = ( y > fitRangeY[0] && 
                                  y < fitRangeY[1]
                                  );
            if (keptNeighbor){
                roo_Meta->setVal(y);
                roo_Weight->setVal(weight);
                // roo_Weight will get overwritten here when adding to RooDataSet but actually does not pick up the value. So we cannot use 
                // roo_Weight.getVal() but the dataset will be weighted: https://root-forum.cern.ch/t/fit-to-a-weighted-unbinned-data-set/33495
                rooData->add(RooArgSet(*roo_Meta,*roo_Weight),weight);
            }
            return keptNeighbor;
        }

        RooArgSet* getParameters(){ 
            // Grabs the parameters of the fig function. Your discriminating variables are not parameters
            return rooSigPlusBkg->getParameters(RooArgList(*roo_Meta));
        }

        void drawFitPlots(vector<vector<float>> &vect, int neighborIdx, TCanvas* c){
            //draw2DPlots(roo_Mpi0, roo_Meta, vect, neighborIdx, eff_nentries, 
            //                rooSigPlusBkg, rooBkg, rooSig, rooData, nsig, nbkg, c);//, model_hist, model_sig, model_bkg);
            draw1DPlots(roo_Meta, vect, neighborIdx, eff_nentries, 
                            rooSigPlusBkg, rooBkg, rooSig, rooData, nsig, nbkg, c);//, model_hist, model_sig, model_bkg);
        }

        RooFitResult* fit(){
            // Looked in multiple places and Moneta always says you cant use Minos with weighted data but should just use SumW2Errors. If I dont
            //   use Minos then nsig+nbkg doesnt even sum to the orignal entries...
            //   lets not use minos for now since it significantly speeds things up. Bootstrapping might help alleviate this problem
            //   https://root-forum.cern.ch/t/parameter-uncertainties-by-rooabspdf-fitto-for-weighted-data-depend-on-minos-option-if-sumw2error-is-set-too/39599/4
            //   https://root.cern.ch/doc/master/rf611__weightedfits_8C.html
            return  rooSigPlusBkg->fitTo(*rooData, 
                    RooFit::Save(), 
                    RooFit::PrintLevel(-1), 
                    RooFit::BatchMode(true), 
                    RooFit::SumW2Error(true)
                    //AsymptoticError(true)), 
                    //Hesse(kFALSE));
                    );
        }
};

#endif
