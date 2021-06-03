#ifndef DRAWPLOTS_H
#define DRAWPLOTS_H

#include "RooAbsRealLValue.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooPolynomial.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TLine.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"

// This is the script that is used to draw the fit results including its x-y projections
// and what the PDFs look like in 2D. 
void draw2DPlots(
        RooAbsRealLValue* x, RooAbsRealLValue* y, vector<vector<float>> vect, int neighborIdx, int nentries, 
        RooAbsPdf* model, RooAbsPdf* bkg, RooAbsPdf* sig, 
        RooDataSet* data, RooAbsRealLValue* nsig, RooAbsRealLValue* nbkg,
        TCanvas* c
        //, TH1* model_hist, TH1* model_sig, TH1* model_bkg
        ){
    using namespace RooFit;
    gStyle->SetLabelSize(0.0525,"xyz");//55,"xyz"); // size of axis value font
    gStyle->SetTitleSize(0.0525,"xyz");//06,"xyz"); // size of axis title font
    gStyle->SetTitleFont(42,"xyz"); // font option
    gStyle->SetPadGridX(1);  // Beni likes this
    gStyle->SetPadGridY(1); 	

    c->Divide(3,2);

    float valX=vect[0][neighborIdx];
    float valY=vect[1][neighborIdx];
    x->setVal(valX);
    y->setVal(valY);
    
    string namex = x->GetName();
    string namey = y->GetName();
    //cout << "x variable name: " << namex << endl;
    //cout << "y variable name: " << namey << endl;

    RooPlot *xframe = x->frame(Title(namex.c_str()));
    RooPlot *xframe2 = y->frame(Title(namey.c_str()));
    TLine* newLine = new TLine();
    newLine->SetLineStyle(2);
    newLine->SetLineWidth(2);

    // Construct plot frame in 'x'
    // Make a second plot frame in x and draw both the
    // data and the pdf in the frame
    data->plotOn(xframe,DataError(RooAbsData::SumW2));
    data->plotOn(xframe2,DataError(RooAbsData::SumW2));

    // Grab binning information
    float locxmin=x->getMin();
    float locxmax=x->getMax();
    float locymin=y->getMin();
    float locymax=y->getMax();
    float locxbins=x->getBins();
    float locybins=y->getBins();
    float locxrange=locxmax-locxmin;
    float locyrange=locymax-locymin;
    float locxbinsize=locxrange/locxbins;
    float locybinsize=locyrange/locybins;
    
    float sigFrac = nsig->getVal()/(nsig->getVal()+nbkg->getVal());

    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////

    // THESE TWO HISTOGRAMS ARE BINNED SO WHEN COMPARING THE PDF VALUES TO THE 
    // HISTOGRAM VALUES WE NEED TO SCALE BY THE BINSIZE

    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    
    // PLOT DATA + PDF INTEGRATED OF X
    float fx;
    float fx_sig;
    float fx_bkg;
    float chiSq;
    RooArgSet* floatPars;
    int nParams;

    c->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetTopMargin(0.225);
    //xframe->GetYaxis()->SetTitleOffset(1);
    model->plotOn(xframe, LineColor(kBlue));
    model->plotOn(xframe, LineColor(kOrange),Components(*bkg));
    model->plotOn(xframe, LineColor(kMagenta),Components(*sig));
    model->paramOn(xframe, Layout(1,1,1)); // essentially hiding the statbox...?
    //model->paramOn(xframe, Layout(0.725,1.0,0.825)); // a decent size if we want to see the statbox
    xframe->getAttText()->SetTextSize(0.01); 
    floatPars = (RooArgSet*)model->getParameters(data)->selectByAttrib("Constant",kFALSE);
    nParams = floatPars->getSize();
    // the order of the objects xframe contains is the order we use plotOn. We can see this if we do xframe->Print()
    chiSq = xframe->chiSquare(xframe->nameOf(1),xframe->nameOf(0),nParams);
    xframe->getAttText()->SetTextSize(0.02); 
    xframe->Draw();
    newLine->SetLineColor(kRed);
    newLine->DrawLine(valX,0,valX,xframe->GetMaximum());    
    fx =  model->getVal(RooArgList(*x))*nentries*locxbinsize;
    newLine->SetLineColor(kBlue);
    newLine->DrawLine(xframe->GetXaxis()->GetXmin(),fx,xframe->GetXaxis()->GetXmax(),fx);
    fx_sig =  sigFrac*sig->getVal(RooArgList(*x))*nentries*locxbinsize;
    newLine->SetLineColor(kMagenta);
    newLine->DrawLine(xframe->GetXaxis()->GetXmin(),fx_sig,xframe->GetXaxis()->GetXmax(),fx_sig);
    fx_bkg =  (1-sigFrac)*bkg->getVal(RooArgList(*x))*nentries*locxbinsize;
    newLine->SetLineColor(kOrange);
    newLine->DrawLine(xframe->GetXaxis()->GetXmin(),fx_bkg,xframe->GetXaxis()->GetXmax(),fx_bkg);
    xframe->SetTitle(("#splitline{Q("+namex+")="+std::to_string(fx_sig/fx)+"}{"+
                        "#splitline{"+"f("+namex+")="+std::to_string(fx)+"   s("+namex+")="+std::to_string(fx_sig)+"}"+
                        "{effective yield="+std::to_string(nentries)+"   reduced ChiSq="+std::to_string(chiSq)+"}}"
                        ).c_str());
    xframe->GetXaxis()->SetTitle((namex+" "+xframe->GetXaxis()->GetTitle()).c_str());

    // PLOT DATA + PDF INTEGRATED OF Y
    c->cd(2);
    gPad->SetTopMargin(0.225);
    xframe2->GetYaxis()->SetTitleOffset(1.6);
    model->plotOn(xframe2, LineColor(kBlue));
    model->plotOn(xframe2, LineColor(kOrange),Components(*bkg));
    model->plotOn(xframe2, LineColor(kMagenta),Components(*sig));
    chiSq = xframe2->chiSquare(xframe2->nameOf(1),xframe2->nameOf(0),nParams);
    xframe2->Draw();
    newLine->SetLineColor(kRed);
    newLine->DrawLine(valY,0,valY,xframe->GetMaximum());    
    fx =  model->getVal(RooArgList(*y))*nentries*locybinsize;
    newLine->SetLineColor(kBlue);
    newLine->DrawLine(xframe2->GetXaxis()->GetXmin(),fx,xframe2->GetXaxis()->GetXmax(),fx);
    fx_sig =  sigFrac*sig->getVal(RooArgList(*y))*nentries*locybinsize;
    newLine->SetLineColor(kMagenta);
    newLine->DrawLine(xframe2->GetXaxis()->GetXmin(),fx_sig,xframe2->GetXaxis()->GetXmax(),fx_sig);
    fx_bkg =  (1-sigFrac)*bkg->getVal(RooArgList(*y))*nentries*locybinsize;
    newLine->SetLineColor(kOrange);
    newLine->DrawLine(xframe2->GetXaxis()->GetXmin(),fx_bkg,xframe2->GetXaxis()->GetXmax(),fx_bkg);
    xframe2->SetTitle(("#splitline{Q("+namey+")="+std::to_string(fx_sig/fx)+"}{f("+namey+")="+std::to_string(fx)+"   s("+namey+")="+std::to_string(fx_sig)+"}").c_str());
    xframe2->SetTitle(("#splitline{Q("+namey+")="+std::to_string(fx_sig/fx)+"}{"+
                        "#splitline{"+"f("+namey+")="+std::to_string(fx)+"   s("+namey+")="+std::to_string(fx_sig)+"}"+
                        "{reduced ChiSq="+std::to_string(chiSq)+"}}"
                        ).c_str());
    xframe2->GetXaxis()->SetTitle((namey+" "+xframe2->GetXaxis()->GetTitle()).c_str());
    
//    delete gROOT->FindObject("");

    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    
    // IN CONTRAST TO THE ABOVE
    // THE FOLLOWING 2D HISTOGRAMS ARE HISTOGRAMS MADE FROM EVALUATING
    // THE PDF ON A GRID. TOTAL = NSIG*SIGNAL + NBKG*BKG WHERE SIGNAL
    // AND BKG ARE BOTH NORMALIZED ALREADY WHEREAS TOTAL IS NORMALIZED
    // TO THE TOTAL ENTRIES.

    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////


    // PLOT PDF ON 2D 
    c->cd(3);
    gPad->SetTopMargin(0.225);
    gPad->SetRightMargin(0.2);
    int pdfBins=locxbins*locybins;
    TH1* model_hist = model->createHistogram((namex+","+namey).c_str(),locxbins,locybins);
    model_hist->SetStats(0);
    model_hist->Scale(locxbinsize*locybinsize);
    newLine->SetLineColor(kRed);
    model_hist->Draw("COLZ");
    newLine->DrawLine(model_hist->GetXaxis()->GetXmin(),valY,model_hist->GetXaxis()->GetXmax(),valY);    
    newLine->DrawLine(valX,model_hist->GetYaxis()->GetXmin(),valX,model_hist->GetYaxis()->GetXmax());    
    fx = model->getVal(RooArgList(*x,*y))*nentries*locxbinsize*locybinsize;
    fx_sig = sig->getVal(RooArgList(*x,*y))*sigFrac*nentries*locxbinsize*locybinsize;
    fx_bkg = bkg->getVal(RooArgList(*x,*y))*(1-sigFrac)*nentries*locxbinsize*locybinsize;
    model_hist->SetTitle(("#splitline{PDF f("+namex+","+namey+")="+std::to_string(fx)+"}{Q("+namex+","+namey+")="+std::to_string(fx_sig/(fx_bkg+fx_sig))+"}").c_str());

    // PLOT SIGNAL PDF ON 2D 
    c->cd(4);
    gPad->SetRightMargin(0.2);
    gPad->SetBottomMargin(0.2);
    TH1* model_sig = sig->createHistogram((namex+","+namey).c_str(),locxbins,locybins);
    model_sig->SetStats(0);
    model_sig->Scale(sigFrac*nentries);//*locxbins/locxrange*locybins/locyrange);
    newLine->SetLineColor(kRed);
    model_sig->Draw("COLZ");
    newLine->DrawLine(model_sig->GetXaxis()->GetXmin(),valY,model_sig->GetXaxis()->GetXmax(),valY);    
    newLine->DrawLine(valX,model_sig->GetYaxis()->GetXmin(),valX,model_sig->GetYaxis()->GetXmax());    
    model_sig->SetTitle(("Signal PDF   f("+namex+","+namey+")="+std::to_string(fx_sig)).c_str());

    // PLOT BACKGROUND PDF ON 2D 
    c->cd(5);
    gPad->SetRightMargin(0.2);
    gPad->SetBottomMargin(0.2);
    TH1* model_bkg = bkg->createHistogram((namex+","+namey).c_str(),locxbins,locybins);
    model_bkg->Scale((1-sigFrac)*nentries);//*locxbins/locxrange*locybins/locyrange);
    model_bkg->SetStats(0);
    newLine->SetLineColor(kRed);
    model_bkg->Draw("COLZ");
    newLine->DrawLine(model_bkg->GetXaxis()->GetXmin(),valY,model_bkg->GetXaxis()->GetXmax(),valY);    
    newLine->DrawLine(valX,model_bkg->GetYaxis()->GetXmin(),valX,model_bkg->GetYaxis()->GetXmax());    
    model_bkg->SetTitle(("Bkg PDF   f("+namex+","+namey+")="+std::to_string(fx_bkg)).c_str());
}

void draw1DPlots(
        RooAbsRealLValue* x, vector<vector<float>> vect, int neighborIdx, int nentries, 
        RooAbsPdf* model, RooAbsPdf* bkg, RooAbsPdf* sig, 
        RooDataSet* data, RooAbsRealLValue* nsig, RooAbsRealLValue* nbkg,
        TCanvas* c
        //, TH1* model_hist, TH1* model_sig, TH1* model_bkg
        ){
    using namespace RooFit;
    gStyle->SetLabelSize(0.0525,"xyz");//55,"xyz"); // size of axis value font
    gStyle->SetTitleSize(0.0525,"xyz");//06,"xyz"); // size of axis title font
    gStyle->SetTitleFont(42,"xyz"); // font option
    gStyle->SetPadGridX(1);  // Beni likes this
    gStyle->SetPadGridY(1); 	

    float valX=vect[0][neighborIdx];
    x->setVal(valX);
    
    string namex = x->GetName();
    cout << "x variable name: " << namex << endl;

    RooPlot *xframe = x->frame(Title(namex.c_str()));
    TLine* newLine = new TLine();
    newLine->SetLineStyle(2);
    newLine->SetLineWidth(2);

    // Construct plot frame in 'x'
    // Make a second plot frame in x and draw both the
    // data and the pdf in the frame
    data->plotOn(xframe,DataError(RooAbsData::SumW2));

    // Grab binning information
    float locxmin=x->getMin();
    float locxmax=x->getMax();
    float locxbins=x->getBins();
    float locxrange=locxmax-locxmin;
    float locxbinsize=locxrange/locxbins;
    
    float sigFrac = nsig->getVal()/(nsig->getVal()+nbkg->getVal());

    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////

    // THESE TWO HISTOGRAMS ARE BINNED SO WHEN COMPARING THE PDF VALUES TO THE 
    // HISTOGRAM VALUES WE NEED TO SCALE BY THE BINSIZE

    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    
    // PLOT DATA + PDF INTEGRATED OF X
    float fx;
    float fx_sig;
    float fx_bkg;
    float chiSq;
    RooArgSet* floatPars;
    int nParams;

    //gPad->SetLeftMargin(0.15);
    //gPad->SetTopMargin(0.225);
    //xframe->GetYaxis()->SetTitleOffset(1);
    model->plotOn(xframe, LineColor(kBlue));
    model->plotOn(xframe, LineColor(kOrange),Components(*bkg));
    model->plotOn(xframe, LineColor(kMagenta),Components(*sig));
    model->paramOn(xframe, Layout(1,1,1)); // essentially hiding the statbox...?
    //model->paramOn(xframe, Layout(0.725,1.0,0.825)); // a decent size if we want to see the statbox
    xframe->getAttText()->SetTextSize(0.01); 
    floatPars = (RooArgSet*)model->getParameters(data)->selectByAttrib("Constant",kFALSE);
    nParams = floatPars->getSize();
    // the order of the objects xframe contains is the order we use plotOn. We can see this if we do xframe->Print()
    chiSq = xframe->chiSquare(xframe->nameOf(1),xframe->nameOf(0),nParams);
    xframe->getAttText()->SetTextSize(0.02); 
    xframe->Draw();
    newLine->SetLineColor(kRed);
    newLine->DrawLine(valX,0,valX,xframe->GetMaximum());    
    fx =  model->getVal(RooArgList(*x))*nentries*locxbinsize;
    newLine->SetLineColor(kBlue);
    newLine->DrawLine(xframe->GetXaxis()->GetXmin(),fx,xframe->GetXaxis()->GetXmax(),fx);
    fx_sig =  sigFrac*sig->getVal(RooArgList(*x))*nentries*locxbinsize;
    newLine->SetLineColor(kMagenta);
    newLine->DrawLine(xframe->GetXaxis()->GetXmin(),fx_sig,xframe->GetXaxis()->GetXmax(),fx_sig);
    fx_bkg =  (1-sigFrac)*bkg->getVal(RooArgList(*x))*nentries*locxbinsize;
    newLine->SetLineColor(kOrange);
    newLine->DrawLine(xframe->GetXaxis()->GetXmin(),fx_bkg,xframe->GetXaxis()->GetXmax(),fx_bkg);
    xframe->SetTitle(("#splitline{Q("+namex+")="+std::to_string(fx_sig/fx)+"}{"+
                        "#splitline{"+"f("+namex+")="+std::to_string(fx)+"   s("+namex+")="+std::to_string(fx_sig)+"}"+
                        "{effective yield="+std::to_string(nentries)+"   reduced ChiSq="+std::to_string(chiSq)+"}}"
                        ).c_str());
    xframe->GetXaxis()->SetTitle((namex+" "+xframe->GetXaxis()->GetTitle()).c_str());
}

#endif

