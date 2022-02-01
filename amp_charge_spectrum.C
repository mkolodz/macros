#include<iostream>
#include<fstream>
#define N_SIPMS 64
#define RESCALE1 1
#define RESCALE2 4.096

#define N_BINS 100
#define AMP_LOW 0.0
#define AMP_HIGH 3000.
#define CH_LOW 0.
#define CH_HIGH 30000.
#define BL_LOW -20.
#define BL_HIGH 20.
#define BLS_LOW 0.
#define BLS_HIGH 20.

// #define AMP_LOW -1000.
// #define AMP_HIGH 1000.
// #define CH_LOW -10000.
// #define CH_HIGH 10000.

/*
#define N_BINS 500
#define AMP_LOW -1000.
#define AMP_HIGH 1000.
#define CH_LOW -10000.
#define CH_HIGH 10000.
#define BL_LOW -1000.
#define BL_HIGH 10000.
#define BLS_LOW -1000.
#define BLS_HIGH 1000.*/


void format_hist(TH1 * h) {
    h->SetTitle("");
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.85);
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetMaxDigits(3);
    h->GetZaxis()->SetLabelSize(0.06);
    gPad->SetLeftMargin(0.125);
    gPad->SetTopMargin(0.07);
    gPad->SetRightMargin(0.125);
    gPad->SetBottomMargin(0.125);
}
Bool_t amp_charge_spectrum(){

    char side = 'l';
    double blSig = 10.0;
    double meas_time = 1800.;
    gStyle->SetPalette(57);
    //gStyle->SetOptStat(0);

    TFile *file = new TFile("sifi_results.root", "READ");
    TTree *tree = (TTree*)file->Get("S");
    
    TString selection_blsig = Form("SDDSamples.data.signal_%c.fBL_sigma*10e3>>htemp3(%d, %f, %f)", side, N_BINS,BLS_LOW,BLS_HIGH);
    TString cut_blsig = "SDDSamples.data.module==0";
    TString selection_bl = Form("SDDSamples.data.signal_%c.fBL*10e3>>htemp4(%d, %f, %f)", side,N_BINS,BL_LOW, BL_HIGH);
    TString cut_bl = "SDDSamples.data.module==0";
    
    tree->Draw(selection_blsig, cut_blsig, "");
    TH1F *hblsig = (TH1F*)gROOT->FindObjectAny("htemp3");
    hblsig->GetXaxis()->SetTitle("baseline sigma [mV]");
    hblsig->SetTitle(" ");
    hblsig->GetYaxis()->SetTitle("counts");
    hblsig->GetStdDev();

    tree->Draw(selection_bl, cut_bl, "");
    TH1F *hbl = (TH1F*)gROOT->FindObjectAny("htemp4");
    hbl->GetXaxis()->SetTitle("baseline [mV]");
    hbl->SetTitle(" ");
    hbl->GetYaxis()->SetTitle("counts");
    
    TString selection_Charge = Form("SDDSamples.data.signal_%c.fCharge*10e3>>htemp5(%d, %f, %f)", side, N_BINS,CH_LOW,CH_HIGH);
    //TString cut_Charge = Form("SDDSamples.data.module==0 && SDDSamples.data.signal_%c.fVeto==0 && SDDSamples.data.signal_%c.fBL_sigma<%f", side, side, blSig);
    TString cut_Charge = "SDDSamples.data.module==0";
    TString selection_Amp = Form("SDDSamples.data.signal_%c.fAmp*10e3>>htemp6(%d, %f, %f)", side, N_BINS,AMP_LOW,AMP_HIGH);
    //TString cut_Amp = Form("SDDSamples.data.module==0 && SDDSamples.data.signal_%c.fVeto==0 && SDDSamples.data.signal_%c.fBL_sigma<%f", side, side, blSig);
    TString cut_Amp = "SDDSamples.data.module==0";
    
    //tree->Draw(selection_Charge, cut_Charge, "");
    tree->Draw(selection_Charge, "", "");
    TH1F *hq = (TH1F*)gROOT->FindObjectAny("htemp5");
    hq->GetXaxis()->SetTitle("charge [pVs]");
    hq->GetYaxis()->SetTitle("counts");
    hq->SetTitle(" ");
    
    tree->Draw(selection_Amp, cut_Amp, "");
    TH1F *ha = (TH1F*)gROOT->FindObjectAny("htemp6");
    ha->GetXaxis()->SetTitle("amplitude [mV]");
    ha->GetYaxis()->SetTitle("counts");
    ha->SetTitle(" ");
    
    TCanvas *can_AmpCharge = new TCanvas("can_AmpCharge", "can_AmpCharge", 1200*1.5, 300*1.5);
    can_AmpCharge->SetRightMargin(0.09);
    can_AmpCharge->SetLeftMargin(0.15);
    can_AmpCharge->SetBottomMargin(0.15);
    can_AmpCharge->Divide(4,1);
    can_AmpCharge->cd(1);
        format_hist(hq);
        gPad->SetGrid(1,1);
        //auto legend = new TLegend(0.5,0.5,0.99,0.75);
        //legend->SetTextSize(0.05);
        //legend->AddEntry((TObject*)0, Form("Rate: %.2f [cps]", ha->GetEntries()/meas_time), "");
        hq->Draw();
        //legend->Draw();
    can_AmpCharge->cd(2);
        format_hist(ha);
        gPad->SetGrid(1,1);
        ha->Draw();
    can_AmpCharge->cd(3); 
        format_hist(hbl);
        gPad->SetGrid(1,1);
        hbl->Draw();
    can_AmpCharge->cd(4);
        format_hist(hblsig);
        gPad->SetGrid(1,1);
        hblsig->Draw();
    can_AmpCharge->cd(5);
   
    can_AmpCharge->Write();
    can_AmpCharge->SaveAs("can_AmpCharge.pdf");

    
    return true;
}
