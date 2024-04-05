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

int checkTimeOfRun(TString path = "/scratch1/gccb/data/Jan2023Beam/root/sortedNoBeam/sortedRun00596_single_b1.root")
{
    
    gROOT->Reset();
    gStyle->SetOptStat(111111);

    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(path);
    if (!f) 
        f = new TFile(path,"READ");
    
    TTree* events = (TTree*)f->Get("events");
    
    //Declaration of leaves types
    Float_t         step1;
    Float_t         step2;
    Long64_t        time;
    UInt_t          channelID;
    Float_t         tot;
    Float_t         energy;
    UShort_t        tacID;
    Int_t           xi;
    Int_t           yi;
    Float_t         x;
    Float_t         y;
    Float_t         z;
    Float_t         tqT;
    Float_t         tqE;

    Long64_t time_first = 0;
    Long64_t time_last = 0;
//     Double_t ps_to_ns = 0.001;
    // Set branch addresses.
    events->SetBranchAddress("step1",&step1);
    events->SetBranchAddress("step2",&step2);
    events->SetBranchAddress("time",&time);
    events->SetBranchAddress("channelID",&channelID);
    events->SetBranchAddress("tot",&tot);
    events->SetBranchAddress("energy",&energy);
    events->SetBranchAddress("tacID",&tacID);
    events->SetBranchAddress("xi",&xi);
    events->SetBranchAddress("yi",&yi);
    events->SetBranchAddress("x",&x);
    events->SetBranchAddress("y",&y);
    events->SetBranchAddress("z",&z);
    events->SetBranchAddress("tqT",&tqT);
    events->SetBranchAddress("tqE",&tqE);
    
    Long64_t nentries = events->GetEntries();

    events->GetEntry(0);
    time_first = time;
//     time_first = time_first*ps_to_ns;
    std::cout << "first timestamp:" << time_first << std::endl;

    events->GetEntry(nentries-1); 
    time_last = time;
//     time_last = time_last*ps_to_ns;
    std::cout << "last timestamp:" << time_last << std::endl;
    std::cout << "last-first = time of measurement [ps]: " << time_last-time_first << std::endl;
    std::cout << "which is, in s: "<< (time_last-time_first)/10e12 << std::endl;
return 0;
}
