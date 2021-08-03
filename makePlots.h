#ifndef MAKEPLOTS_H
#define MAKEPLOTS_H

#include <ctime>
#include <math.h> 
#include "configSettings.h"

using namespace std;

// This will stack our histograms and make them pretty. 
void makeStackedHist(TH1F* truth, TH1F* tot, TH1F* sig, TH1F* bkg, TH1F* sig_sb, TH1F* bkg_sb, string name,string baseDir){
    	TCanvas *allCanvases = new TCanvas(name.c_str(),"",1440,900);
        allCanvases->Divide(2,2);
	TLegend* leg1 = new TLegend(0.8,0.8,1.0,1.0);
	TLegend* leg1_sb = new TLegend(0.8,0.8,1.0,1.0);
	TLegend* leg1_overlay = new TLegend(0.8,0.8,1.0,1.0);
	THStack* stackedHists = new THStack("stackedHists","");
	THStack* stackedHists_sb = new THStack("stackedHists_sb","");
	THStack* stackedHists_overlay = new THStack("stackedHists_overlay","");

        bool haveTruth=truth->GetEntries()>0;

	bkg->SetFillColorAlpha(kViolet-5,0.3);
	bkg->SetLineColorAlpha(kViolet-5,0);
	sig->SetLineColorAlpha(kRed+1,1);
	sig->SetLineWidth(2);

	bkg_sb->SetFillColorAlpha(kViolet-5,0.3);
	bkg_sb->SetLineColorAlpha(kViolet-5,0);
	sig_sb->SetLineColorAlpha(kBlue+1,1);
	sig_sb->SetLineWidth(2);

	leg1->AddEntry(bkg,"Q_Bkg","f");
	leg1->AddEntry(sig,"Q_Sig","l");
	leg1->AddEntry(tot,"Tot","l");
	leg1_sb->AddEntry(bkg_sb,"SB_Bkg","f");
	leg1_sb->AddEntry(sig_sb,"SB_Sig","l");
	leg1_sb->AddEntry(tot,"Tot","l");
	leg1_overlay->AddEntry(sig,"Q_Sig","l");
	leg1_overlay->AddEntry(sig_sb,"SB_Sig","l");

        if (haveTruth){
            truth->SetFillColorAlpha(kGray,0.7);
            truth->SetLineColorAlpha(8,0);
            leg1_overlay->AddEntry(truth,"matchedThrown","f");
        }

        allCanvases->cd(1);
	stackedHists->Add(tot,"HIST");
	stackedHists->Add(bkg,"HIST");
	stackedHists->Add(sig,"HIST");
	stackedHists->Draw("nostack");
        stackedHists->SetMaximum(round(1.05*tot->GetMaximum()));
	leg1->Draw();
	stackedHists->GetXaxis()->SetTitle(tot->GetXaxis()->GetTitle());
	stackedHists->GetYaxis()->SetTitle(tot->GetYaxis()->GetTitle());

        allCanvases->cd(2);
	stackedHists_sb->Add(tot,"HIST");
	stackedHists_sb->Add(bkg_sb,"HIST");
	stackedHists_sb->Add(sig_sb,"HIST");
	stackedHists_sb->Draw("nostack");
        stackedHists_sb->SetMaximum(round(1.05*tot->GetMaximum()));
	leg1_sb->Draw();
	stackedHists_sb->GetXaxis()->SetTitle(tot->GetXaxis()->GetTitle());
	stackedHists_sb->GetYaxis()->SetTitle(tot->GetYaxis()->GetTitle());

        allCanvases->cd(3);
        if (haveTruth)
            stackedHists_overlay->Add(truth,"HIST");
        stackedHists_overlay->Add(sig,"HIST");
        stackedHists_overlay->Add(sig_sb,"HIST");
	stackedHists_overlay->Draw("nostack");
        stackedHists_overlay->SetMaximum(round(1.05*sig->GetMaximum()));
	stackedHists_overlay->GetXaxis()->SetTitle(tot->GetXaxis()->GetTitle());
	stackedHists_overlay->GetYaxis()->SetTitle(tot->GetYaxis()->GetTitle());
        leg1_overlay->Draw();

	allCanvases->SaveAs((baseDir+"/"+name+".png").c_str());
}

void make2DHistsOnPads(TH2F* truth, TH2F* tot, TH2F* sig, TH2F* bkg, TH2F* sig_sb, TH2F* bkg_sb, string name,string baseDir){
    	TCanvas *allCanvases = new TCanvas(name.c_str(),"",1440,900);
        allCanvases->Divide(3,2);
	TLegend* leg1 = new TLegend(0.8,0.8,1.0,1.0);
	TLegend* leg1_sb = new TLegend(0.8,0.8,1.0,1.0);
	TLegend* leg1_overlay = new TLegend(0.8,0.8,1.0,1.0);

        allCanvases->cd(1);
        tot->Draw("COLZ");

        allCanvases->cd(2);
        sig->Draw("COLZ");

        allCanvases->cd(3);
        bkg->Draw("COLZ");

        allCanvases->cd(5);
        sig_sb->Draw("COLZ");

        allCanvases->cd(6);
        bkg_sb->Draw("COLZ");

	allCanvases->SaveAs((baseDir+"/"+name+".png").c_str());
}

#endif


