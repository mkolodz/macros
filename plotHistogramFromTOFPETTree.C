#include "TH1D.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include <RtypesCore.h>

#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <ctime>
#include <algorithm>

int plotHistogramFromTOFPETTree(TString path1 = "/scratch1/gccb/data/Jan2023Beam/root/sortedNoBeam/sortedRun00596_single_b0.root",TString path2 = "/scratch1/gccb/data/Jan2023Beam/root/run00596_single_nlc.root", TString path3 = "/scratch1/gccb/data/Jan2023Beam/root/run00596_single.root")
{
    
    gROOT->Reset();
    gStyle->SetOptStat(111111);

    TFile *f1 = (TFile*)gROOT->GetListOfFiles()->FindObject(path1);
    if (!f1) f1 = new TFile(path1,"READ");
    TTree* events1 = (TTree*)f1->Get("events");
    
    TFile *f2 = (TFile*)gROOT->GetListOfFiles()->FindObject(path2);
    if (!f2) f2 = new TFile(path2, "READ");
    TTree* events2 = (TTree*)f2->Get("events");
    
    TFile *f3 = (TFile*)gROOT->GetListOfFiles()->FindObject(path3);
    if (!f3) f3 = new TFile(path3,"READ");
    TTree* events3 = (TTree*)f3->Get("events");
    
    TH1D * hEnergy1 = new TH1D("hEnergy1", "hEnergy1", 100, -30,150); 
    TH1D * hEnergy2 = new TH1D("hEnergy2", "hEnergy2", 100, -30,150);  
    TH1D * hEnergy3 = new TH1D("hEnergy3", "hEnergy3", 100, -30,150);  
//     Long64_t        time;
//     UInt_t          channelID;
    Float_t         energy1, energy2, energy3;

//     events1->SetBranchAddress("time",&time);
//     events1->SetBranchAddress("channelID",&channelID);
    events1->SetBranchAddress("energy",&energy1);
    events2->SetBranchAddress("energy",&energy2);
    events3->SetBranchAddress("energy",&energy3);
    Long64_t nentries1 = events1->GetEntries();
    Long64_t nentries2 = events2->GetEntries();
    Long64_t nentries3 = events3->GetEntries();
//     for(int i=0; i<nentries1; i++){
    for(int i=0; i<20000; i++){
        events1->GetEntry(i);
        hEnergy1->Fill(energy1);
//         std::cout << energy1 << std::endl;
    }
//     std::cout << std::endl;
//         for(int i=0; i<nentries2; i++){
    for(int i=0; i<20000; i++){
        events2->GetEntry(i);
        hEnergy2->Fill(energy2);
//         std::cout << energy2 << std::endl;
    }
//     std::cout << std::endl;
//         for(int i=0; i<nentries3; i++){
    for(int i=0; i<20000; i++){
        events3->GetEntry(i);
        hEnergy3->Fill(energy3);
//         std::cout << energy3 << std::endl;
    }
//     std::cout << std::endl;
//     std::cout << std::endl;
    TCanvas * c = new TCanvas("c", "c", 600, 600);
    c->cd();
    hEnergy1->Scale(1./20000);
    hEnergy1->Draw();
    hEnergy1->SetLineColor(kRed);
    hEnergy2->Scale(1./20000);
    hEnergy2->Draw("SAME");
    hEnergy2->SetLineColor(kBlue);
    hEnergy3->Scale(1./20000);
    hEnergy3->Draw("SAME");
    hEnergy3->SetLineColor(kGreen);
    gPad->BuildLegend();
    
    
    
//     time_first = time;
// //     time_first = time_first*ps_to_ns;
//     std::cout << "first timestamp:" << time_first << std::endl;
// 
//     events1->GetEntry(nentries1-1); 
//     time_last = time;
// //     time_last = time_last*ps_to_ns;
//     std::cout << "last timestamp:" << time_last << std::endl;
//     std::cout << "last-first = time of measurement [ps]: " << time_last-time_first << std::endl;
//     std::cout << "which is, in s: "<< (time_last-time_first)/10e12 << std::endl;
return 0;
}
