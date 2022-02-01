#include<iostream>
#include<fstream>
#define N_SIPMS 64
#define RESCALE1 1
#define RESCALE2 4.096
#define N_BINS 100
#define AMP_LOW 20.0
#define AMP_HIGH 100.
#define CH_LOW 600.
#define CH_HIGH 2500.
#define BL_LOW -50.
#define BL_HIGH 50.
#define BLS_LOW 0.0
#define BLS_HIGH 100.

// #define N_BINS 100
// #define AMP_LOW 0.0
// #define AMP_HIGH 3000.
// #define CH_LOW 0.
// #define CH_HIGH 30000.
// #define BL_LOW -20.
// #define BL_HIGH 20.
// #define BLS_LOW 0.
// #define BLS_HIGH 20.

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
Bool_t amp_charge_spectrum_automatic_range(){

    char side = 'l';
    double blSig = 10.0;
    double meas_time = 1800.;
    gStyle->SetPalette(57);
    gStyle->SetOptStat(0);

    TFile *file = new TFile("/scratch/gccb/data/202110/scope/2021_10_27/sifi_results_35_42_merged.root", "READ");
    TTree *tree = (TTree*)file->Get("S");
    TFile *file35 = new TFile("/scratch/gccb/data/202110/scope/2021_10_27/sifi_results_35_amp.root", "READ");
    TTree *tree35 = (TTree*)file35->Get("S");
    TFile *file42 = new TFile("/scratch/gccb/data/202110/scope/2021_10_27/sifi_results_42_amp.root", "READ");
    TTree *tree42 = (TTree*)file42->Get("S");
    
    
    TString selection = Form("SDDSamples.data.signal_%c.fAmp*1000:SDDSamples.data.signal_%c.fCharge*1000>>htemp1(%d,%f, %f, %d, %f, %f)", side, side, N_BINS, CH_LOW,CH_HIGH,N_BINS,AMP_LOW,AMP_HIGH);
    TString selection_all = Form("SDDSamples.data.signal_%c.fAmp*1000:SDDSamples.data.signal_%c.fCharge*1000>>htemp2(%d,%f, %f, %d, %f, %f)", side, side, N_BINS, CH_LOW,CH_HIGH,N_BINS,AMP_LOW,AMP_HIGH);
    TString cut = Form("SDDSamples.data.module==0 && SDDSamples.data.signal_%c.fVeto==0 && SDDSamples.data.signal_%c.fBL_sigma<%f", side, side, blSig);
//         TString cut = Form("SDDSamples.data.module==0 && SDDSamples.data.signal_%c.fBL_sigma<%f", side, blSig);
    TString cut_all = "SDDSamples.data.module==0";
    
    tree->Draw(selection, cut, "");
    TH2F *hcut = (TH2F*)gROOT->FindObjectAny("htemp1");
    hcut->GetXaxis()->SetTitle("charge [pVs]");
    hcut->GetYaxis()->SetTitle("amplitude [mV]");
    hcut->SetTitle(" ");
    
    tree->Draw(selection_all, cut_all, "");
    TH2F *hall = (TH2F*)gROOT->FindObjectAny("htemp2");
    //hall->SetCanExtend(TH2::kAllAxes);
    hall->GetXaxis()->SetTitle("charge [pVs]");
    hall->GetYaxis()->SetTitle("amplitude [mV]");
    hall->SetTitle(" ");

    TString selection_Charge = Form("SDDSamples.data.signal_%c.fCharge*1000>>htemp5(%d, %f, %f)", side, N_BINS,CH_LOW,CH_HIGH);
    //TString cut_Charge = Form("SDDSamples.data.module==0 && SDDSamples.data.signal_%c.fVeto==0 && SDDSamples.data.signal_%c.fBL_sigma<%f", side, side, blSig);
    TString cut_Charge = "SDDSamples.data.module==0";
    TString selection_Amp = Form("SDDSamples.data.signal_%c.fAmp*1000>>htemp6(%d, %f, %f)", side, N_BINS,AMP_LOW,AMP_HIGH);
    //TString cut_Amp = Form("SDDSamples.data.module==0 && SDDSamples.data.signal_%c.fVeto==0 && SDDSamples.data.signal_%c.fBL_sigma<%f", side, side, blSig);
    TString cut_Amp = "SDDSamples.data.module==0";
 
    TString selection_Charge_35 = Form("SDDSamples.data.signal_%c.fCharge*1000>>htemp7(%d, %f, %f)", side, N_BINS,CH_LOW,CH_HIGH);
    TString selection_Charge_42 = Form("SDDSamples.data.signal_%c.fCharge*1000>>htemp8(%d, %f, %f)", side, N_BINS,CH_LOW,CH_HIGH);
    TString selection_Amp_35 = Form("SDDSamples.data.signal_%c.fAmp*1000>>htemp9(%d, %f, %f)", side, N_BINS,AMP_LOW,AMP_HIGH);
    TString selection_Amp_42 = Form("SDDSamples.data.signal_%c.fAmp*1000>>htemp10(%d, %f, %f)", side, N_BINS,AMP_LOW,AMP_HIGH);
    
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
    
    tree35->Draw(selection_Charge_35, "", "");
    TH1F *hq35 = (TH1F*)gROOT->FindObjectAny("htemp7");
    tree42->Draw(selection_Charge_42, "", "");
    TH1F *hq42 = (TH1F*)gROOT->FindObjectAny("htemp8");
    tree35->Draw(selection_Amp_35, "", "");
    TH1F *ha35 = (TH1F*)gROOT->FindObjectAny("htemp9");
    tree42->Draw(selection_Amp_42, "", "");
    TH1F *ha42 = (TH1F*)gROOT->FindObjectAny("htemp10");    
    
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
//     
//     TString selection_Charge = Form("SDDSamples.data.signal_%c.fCharge*10e3>>htemp5(%d, %f, %f)", side, N_BINS,CH_LOW,CH_HIGH);
//     TString cut_Charge = Form("SDDSamples.data.module==0 && SDDSamples.data.signal_%c.fVeto==0 && SDDSamples.data.signal_%c.fBL_sigma<%f", side, side, blSig);
//     TString cut_Charge = "SDDSamples.data.module==0";
//     TString selection_Amp = Form("SDDSamples.data.signal_%c.fAmp*10e3>>htemp6(%d, %f, %f)", side, N_BINS,AMP_LOW,AMP_HIGH);
//     TString cut_Amp = Form("SDDSamples.data.module==0 && SDDSamples.data.signal_%c.fVeto==0 && SDDSamples.data.signal_%c.fBL_sigma<%f", side, side, blSig);
//     TString cut_Amp = "SDDSamples.data.module==0";
    
    //tree->Draw(selection_Charge, cut_Charge, "");
    //tree->Draw(selection_Charge, "", "");
    //tree->Draw("SDDSamples.data.signal_l.fCharge");
    //TH1F *hq = (TH1F*)gROOT->FindObjectAny("htemp5");
    //hq->GetXaxis()->SetTitle("charge [pVs]");
    //hq->GetYaxis()->SetTitle("counts");
    //hq->SetTitle(" ");
/*    
    tree->Draw(selection_Amp, cut_Amp, "");
    TH1F *ha = (TH1F*)gROOT->FindObjectAny("htemp6");
    ha->GetXaxis()->SetTitle("amplitude [mV]");
    ha->GetYaxis()->SetTitle("counts");
    ha->SetTitle(" ");*/
    
    TCanvas *can_AmpCharge = new TCanvas("can_AmpCharge", "can_AmpCharge", 800*1.5, 1200*1.5);
    can_AmpCharge->SetRightMargin(0.09);
    can_AmpCharge->SetLeftMargin(0.15);
    can_AmpCharge->SetBottomMargin(0.15);
    can_AmpCharge->Divide(2,3);
    can_AmpCharge->cd(1);
    hq->Draw();
    hq35->SetLineColor(kRed);
    hq35->Draw("same");    
    hq42->SetLineColor(kGreen);
    hq42->Draw("same");
        //format_hist(hq);
        gPad->SetGrid(1,1);
        //auto legend = new TLegend(0.5,0.5,0.99,0.75);
        //legend->SetTextSize(0.05);
        //legend->AddEntry((TObject*)0, Form("Rate: %.2f [cps]", ha->GetEntries()/meas_time), "");
        //hq->Draw();
        //legend->Draw();
        //tree->Draw("SDDSamples.data.signal_l.fCharge");

    auto legend = new TLegend(0.6,0.8,0.99,0.99); // option "C" allows to center the header
    legend->AddEntry(hq35,"THR=35mV","l"); //p stands for "points"
    legend->AddEntry(hq42,"THR=42mV","l");
    legend->AddEntry(hq,"merged data with both THRs","l");
    legend->Draw();
    can_AmpCharge->cd(2);
    ha->Draw();
    ha35->SetLineColor(kRed);
    ha35->Draw("same");
    ha42->SetLineColor(kGreen);
    ha42->Draw("same");
    auto legend2 = new TLegend(0.6,0.8,0.99,0.99); // option "C" allows to center the header
    legend2->AddEntry(ha35,"THR=35mV","l"); //p stands for "points"
    legend2->AddEntry(ha42,"THR=42mV","l");
    legend2->AddEntry(ha,"merged data with both THRs","l");
    legend2->Draw();
    //tree->Draw("SDDSamples.data.signal_l.fAmp");
//         format_hist(ha);
        gPad->SetGrid(1,1);
//         ha->Draw();
    can_AmpCharge->cd(3); 
        format_hist(hall);
        gPad->SetGrid(1,1);
        hall->Draw("colz");

    
    can_AmpCharge->cd(4);
        tree->Draw("SDDSamples.data.signal_l.fBL");
        format_hist(hbl);
        gPad->SetGrid(1,1);
        hbl->Draw();
    
    can_AmpCharge->cd(5);    
        tree->Draw("SDDSamples.data.signal_l.fBL_sigma");
        format_hist(hblsig);
        gPad->SetGrid(1,1);
        hblsig->Draw();

    can_AmpCharge->cd(6);
   Int_t n = 4;
   Double_t x[4] = {1.,2.,3.,4.};
    Double_t y[4] = {1.02948e+03,1.20311e+03,  1.38952e+03, 1.56912e+03};
   TGraph* gr = new TGraph(n,x,y);
//    Double_t fun(){
//        return 
// }
   gr->Fit("pol1");

   //std::cout <<    gr->GetParameter(0) << std::endl;
        //gr->SetMaximum(1000000.); 
        //gr->SetMinimum(1.);
        gr->Draw("AP*");
        gr->SetTitle(";#PE (arbitrary) ;Charge [pVs]");
        //gr->SetMarkerColor(4);

    
    can_AmpCharge->Write();
    can_AmpCharge->SaveAs("can_AmpCharge.pdf");

    
    return true;
}
