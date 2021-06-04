#ifndef GETSBWEIGHT_H
#define GETSBWEIGHT_H

#include <iostream>
#include "./auxilliary/customPDFs/bivariateGaus/bivariateGaus.h"
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
#include <RooClassFactory.h>
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
        float initMassX = 0.135059;
        float initMassY = 0.549354;
        float initSigmaX = 0.00612767;
        float initSigmaY = 0.0166494;
        float initRho = 0.0;
        float initBernA = 0.5;
        float initBernB = 0.5;
        float initBernC = 0.5;
        float initBernD = 0.5;
        float initBernE = 0.5;
        float initbkgPeakFrac=0.5;
        std::vector<float> fitRangeY={0.36,0.75};
        std::vector<float> fitRangeX={0.085,0.185};

        ///////////////////////////////
        // DEFINE VARIABLES FOR DATASET
        ///////////////////////////////
        float eff_nentries;
        RooRealVar* y;
        RooRealVar* x;
        RooRealVar* w;
        RooDataSet* rooData;
        ////////////////////////
        // DEFINE THE SIGNAL PDF
        ////////////////////////
        RooRealVar* px;
        RooRealVar* sx;
        RooRealVar* py;
        RooRealVar* sy;
        RooRealVar* rho;
        RooGaussian* rooGausPi0;
        RooGaussian* rooGausEta;
        //RooGenericPdf* rooSig;
        //RooProdPdf* rooSig;
        //bivariateGaus* rooSig;
        RooAbsPdf* rooSig;
        ////////////////////////
        // DEFINE THE BKG PDF
        ////////////////////////
        RooRealVar* bern_parA;
        RooRealVar* bern_parB;
        RooRealVar* bern_parC;
        RooRealVar* bern_parD;
        RooGenericPdf* rooBkgX;
        RooGenericPdf* rooBkgY;
        RooGaussian* rooGausPi0_bkg;
        RooRealVar* bkgPeakFrac;
        RooAddPdf* rooBkgXplusPi0Peak;
        RooProdPdf* rooBkg;
        ////////////////////////
        // DEFINE THE TOTAL PDF
        ////////////////////////
        RooRealVar* nsig;
        RooRealVar* nbkg;
        RooAddPdf* rooSigBkg;
        
	fitManager(string iProcess){ 
            /////////////// VARIABLES
            y = new RooRealVar{("y"+iProcess).c_str(),"Mass GeV",(fitRangeY[1]-fitRangeY[0])/2+fitRangeY[0],fitRangeY[0],fitRangeY[1]};
            y->setRange(("roo_fitRangeMeta_"+iProcess).c_str(),fitRangeY[0], fitRangeY[1]);
            y->setBins(100);
            x = new RooRealVar{("x"+iProcess).c_str(),"Mass GeV",(fitRangeX[1]-fitRangeX[0])/2+fitRangeX[0],fitRangeX[0],fitRangeX[1]};
            x->setRange(("roo_fitRangeMpi0_"+iProcess).c_str(),fitRangeX[0], fitRangeX[1]);
            x->setBins(100);
            w = new RooRealVar{("w"+iProcess).c_str(), "Weight", 0, -10, 10}; // Weights can take a wide range
            rooData = new RooDataSet{("rooData"+iProcess).c_str(),"rooData",RooArgSet(*x,*y,*w),RooFit::WeightVar(*w)};
            /////////////// FOR SIGNAL PDF
            px = new RooRealVar{("px"+iProcess).c_str(),"px",initMassX};
            sx = new RooRealVar{("sx"+iProcess).c_str(),"sx",initSigmaX,initSigmaX,initSigmaX*3};
            py = new RooRealVar{("py"+iProcess).c_str(),"py",initMassY};
            sy = new RooRealVar{("sy"+iProcess).c_str(),"sy",initSigmaY,initSigmaY,initSigmaY*3};
            rho = new RooRealVar{("rho"+iProcess).c_str(),"rho",initRho,-0.9,0.9};
            //rooGausPi0 = new RooGaussian{("rooGausPi0_"+iProcess).c_str(), "rooGausPi0", *x, *px, *sx};
            //rooGausEta = new RooGaussian{("rooGausEta_"+iProcess).c_str(), "rooGausEta", *y, *py, *sy};
            /////////////// FOR BKG PDF
            bern_parA = new RooRealVar{("bern_parA"+iProcess).c_str(),"bern_parA",initBernA,0,1};
            bern_parB = new RooRealVar{("bern_parB"+iProcess).c_str(),"bern_parB",initBernB,0,1};
            bern_parC = new RooRealVar{("bern_parC"+iProcess).c_str(),"bern_parC",initBernC,0,1};
            bern_parD = new RooRealVar{("bern_parD"+iProcess).c_str(),"bern_parD",initBernD,0,1};
            rooBkgY = new RooGenericPdf{("rooBkgY"+iProcess).c_str(), "rooBkgY"
                ,("bern_parC"+iProcess+"*y"+iProcess+"+bern_parD"+iProcess+"*(1-y"+iProcess+")").c_str()
                ,RooArgSet(*bern_parC,*bern_parD,*y)};
            rooGausPi0_bkg = new RooGaussian{("rooGausPi0_bkg"+iProcess).c_str(), "rooGausPi0_bkg", *x, *px, *sx};
            rooBkgX = new RooGenericPdf{("rooBkgX"+iProcess).c_str(), "rooBkgX"
                ,("bern_parA"+iProcess+"*x"+iProcess+"+bern_parB"+iProcess+"*(1-x"+iProcess+")").c_str()
                ,RooArgSet(*bern_parA,*bern_parB,*x)};
            bkgPeakFrac = new RooRealVar{("bkgPeakFrac"+iProcess).c_str(),"bkgPeakFrac",initbkgPeakFrac,0,1};
            rooBkgXplusPi0Peak = new RooAddPdf{("rooBkgXplusPi0Peak"+iProcess).c_str(), "rooBkgXplusPi0Peak", RooArgList(*rooGausPi0_bkg,*rooBkgX),RooArgSet(*bkgPeakFrac)};
            /////////////// FOR YIELDS
            nsig = new RooRealVar{("nsig"+iProcess).c_str(),"nsig",(float)kDim/2,0,(float)kDim};
            nbkg = new RooRealVar{("nbkg"+iProcess).c_str(),"nbkg",(float)kDim/2,0,(float)kDim};
            /////////////// CREATING FINAL PDFS
            // -----------------------------------------------------------------
            // HOW TO MAKE YOUR OWN CUSTOM PDFs : MUCH FASTER THAN ROOGENERICPDF
            // 1. Use RooClassFactory to create custom PDF class. Should create a .so file
            // 2. Include the header file here
            // 3. include something like -pathToSoFile when you compile the main program
            // -----------------------------------------------------------------
            //rooSig = new RooProdPdf{("rooSig"+iProcess).c_str(), "rooSig", RooArgSet(*rooGausPi0,*rooGausEta)};
            //rooSig = new RooGenericPdf(("rooSig"+iProcess).c_str(), "Bivariate Gaussian", expression.c_str(),
            //    RooArgList(*x, *y, *rho, *px, *sx, *py, *sy) );
            //string coeff="1/(2*3.14157*sx"+iProcess+"*sy"+iProcess+"*sqrt(1-pow(rho"+iProcess+",2)))";
            //string argCoeff="-1/(2*(1-pow(rho"+iProcess+",2)))";
            //string xterm="((x"+iProcess+"-px"+iProcess+")/sx"+iProcess+")";
            //string yterm="((y"+iProcess+"-py"+iProcess+")/py"+iProcess+")";
            //string crossterm="-2*rho"+iProcess+"*"+xterm+"*"+yterm;
            //string expression=coeff+"*exp("+argCoeff+"*(pow("+xterm+",2)"+crossterm+"+pow("+yterm+",2)))";
            //cout << "Signal expression: " << expression << endl;
            rooSig = new bivariateGaus{("rooSig"+iProcess).c_str(), "Bivariate Gaussian", *x, *y, *px, *sx, *py, *sy, *rho };
            //rooSig =
            //   RooClassFactory::makePdfInstance(("rooSig"+iProcess).c_str(), coeff.c_str(), RooArgList(*x, *y, *rho, *px, *sx, *py, *sy) );
            //rooSig = new RooGenericPdf(("rooSig"+iProcess).c_str(), "Bivariate Gaussian", expression.c_str(),
            //    RooArgList(*x, *y, *rho, *px, *sx, *py, *sy) );
            rooBkg = new RooProdPdf{("rooBkg"+iProcess).c_str(),"rooBkg",RooArgList(*rooBkgXplusPi0Peak,*rooBkgY)};
            rooSigBkg = new RooAddPdf{("rooSumPdf"+iProcess).c_str(), "rooSumPdf", RooArgList(*rooSig,*rooBkg),RooArgSet(*nsig,*nbkg)};
        }

        void reinitialize(float initSigFrac){
            px->setVal(initMassX);
            py->setVal(initMassY);
            sx->setVal(initSigmaX);
            sy->setVal(initSigmaY);
            rho->setVal(initRho);
            bern_parA->setVal(initBernA);
            bern_parB->setVal(initBernB);
            bern_parC->setVal(initBernC);
            bern_parD->setVal(initBernD);
            bkgPeakFrac->setVal(initbkgPeakFrac);
            eff_nentries = rooData->sumEntries();
            nsig->setVal(initSigFrac*eff_nentries);
            nbkg->setVal((1-initSigFrac)*eff_nentries);
        }

        float calculate_q(vector<vector<float>> &vect, int ientry){
            // vect contains the discriminating variables (i.e. Mpi0 and Meta)
            // ientry is the entry in the root tree we are trying to calculate the q-value for
            x->setVal(vect[0][ientry]);
            y->setVal(vect[1][ientry]);
            float sigFrac = nsig->getVal()/(nsig->getVal()+nbkg->getVal());
            float sigPdfVal = sigFrac*rooSig->getVal(RooArgSet(*x,*y));
            float bkgPdfVal = (1-sigFrac)*rooBkg->getVal(RooArgSet(*x,*y));
            float totPdfVal = rooSigBkg->getVal(RooArgSet(*x,*y));
            float qvalue = sigPdfVal/(sigPdfVal+bkgPdfVal);
            //if ((sigPdfVal==0)*(bkgPdfVal==0)*(sigFrac==1)){
            //    qvalue=0;    
            //}
            cout << "\tpostFit(Q=" << qvalue << ") - sigFrac: " << sigFrac << " || sigPdfVal: " << sigPdfVal << " || bkgPdfVal: " << bkgPdfVal << endl;
            return qvalue;
        }

        bool insert(vector<vector<float>> &vect, int ientry, float weight){
            // with 600 neighbors it seems like the fitTo command takes ~2x longer when using Range() argument which selects the fit range.
            // This is equivalent to shrinking the dataset range and fitting over the full range which will save time.
            // It might be useful to think of setting a fit range as to lower the effective number of nearest neighbors. Another thing that lowers
            // the effective number of neighbors is any weights we apply to the filling of the histograms
            float xval = vect[0][ientry];
            float yval = vect[1][ientry];
            bool keptNeighbor = ( yval > fitRangeY[0] && 
                                  yval < fitRangeY[1] &&
                                  xval > fitRangeX[0] &&
                                  xval < fitRangeX[1] );
            if (keptNeighbor){
                x->setVal(xval);
                y->setVal(yval);
                w->setVal(weight);
                // w will get overwritten here when adding to RooDataSet but actually does not pick up the value. So we cannot use 
                // w.getVal() but the dataset will be weighted: https://root-forum.cern.ch/t/fit-to-a-weighted-unbinned-data-set/33495
                rooData->add(RooArgSet(*x,*y,*w),weight);
            }
            return keptNeighbor;
        }

        RooArgSet* getParameters(){ 
            // Grabs the parameters of the fig function. Your discriminating variables are not parameters
            return rooSigBkg->getParameters(RooArgList(*x,*y));
        }

        void drawFitPlots(vector<vector<float>> &vect, int ientry, TH1* dHist_qvaluesBS, double best_qvalue, int iBS, TCanvas* c, TFile* qHistsFile){
            // vect contains the discriminating variables (i.e. Mpi0 and Meta)
            // ientry is the entry in the root tree we are trying to calculate the q-value for
            // dHist_qvaluesBS is the histogram of the boostrapped qvalues
            // best_qvalue is taken from the best fitted qvalue if we do multiple fits per entry
            // iBS is the current bootstrap iteration. Useful if we wanted to see how each bootstrap iteration looks like
            // qHistsFile is the TFile we will save to
            draw2DPlots(x, y, vect, ientry, eff_nentries, 
                            rooSigBkg, rooBkg, rooSig, rooData, nsig, nbkg, 
                            dHist_qvaluesBS, best_qvalue, iBS, c, qHistsFile);
            //draw1DPlots(y, vect, ientry, eff_nentries, 
            //                rooSigBkg, rooBkg, rooSig, rooData, nsig, nbkg, 
            //                dHist_qvaluesBS, best_qvalue, iBS, c, qHistsFile);
        }

        RooFitResult* fit(){
            // Looked in multiple places and Moneta always says you cant use Minos with weighted data but should just use SumW2Errors. If I dont
            //   use Minos then nsig+nbkg doesnt even sum to the orignal entries...
            //   lets not use minos for now since it significantly speeds things up. Bootstrapping might help alleviate this problem
            //   https://root-forum.cern.ch/t/parameter-uncertainties-by-rooabspdf-fitto-for-weighted-data-depend-on-minos-option-if-sumw2error-is-set-too/39599/4
            // SumW2Error - errors calculated from the inverse hessian does not give good coverage for weighted data. Setting to false is better for weighted I think? 
            // AsymptoticError - correct in the large N limit
            //   https://root.cern.ch/doc/master/rf611__weightedfits_8C.html
            return  rooSigBkg->fitTo(*rooData, 
                    RooFit::Save(), 
                    RooFit::PrintLevel(-1), 
                    RooFit::BatchMode(true), 
                    RooFit::SumW2Error(true),//false),
                    //RooFit::AsymptoticError(true), 
                    RooFit::Hesse(kFALSE)
                    );
        }
};

#endif
