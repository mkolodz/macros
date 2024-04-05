// #include "SCategoryManager.h"
// #include "SLoop.h"
// #include "SDDSamples.h"
// #include "SFibersRaw.h"
// #include "SFibersIdentification.h"
// #include "SSiPMHit.h"
// #include "SSiPMCluster.h"
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

int postprocessHeatMaps(TString path = "/scratch1/gccb/data/Jan2023Beam/results/00579_sifi_testing_HEATMAPS_ONLY.root", TString stats = "/scratch1/gccb/data/Jan2023Beam/root/preselectionTimes/run00579_single_nlc_spills.root")
{
    TFile *f = new TFile(path);
    f->ls();
    TH2D * hist = (TH2D*)f->Get("hFiberHeatMap");
    TFile *f2 = new TFile(stats);
    f2->ls();
//     hist->Draw("E");
    TH1D * hist2 = (TH1D*)f2->Get("hstats");
    hist2->Draw();
    path.ReplaceAll(".root", "_wstats.root");
    TFile *output = new TFile(path,"RECREATE");
	output->cd();
    hist->SetDirectory(gDirectory);
    hist->Write();
    hist2->SetDirectory(gDirectory);
    hist2->Write();
//     hist->AddBinContent(1, 3.); 
//     hist->SetEntries(hist->GetEntries());
//     std::cout << hist->GetBinError(10,5) << std::endl;
//     hist->Sumw2();
//     std::cout << hist->GetBinError(1) << std::endl;
    
    return 0;
}
