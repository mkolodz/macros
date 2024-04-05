// #include "/home/magda/project/sifi-framework/sifi-framework-install/include/SCategoryManager.h"
#include "SCategoryManager.h"
#include "SLoop.h"
#include "SDDSamples.h"
#include "SFibersRaw.h"
R__ADD_INCLUDE_PATH(/home/magda/project/sifi-framework/sifi-framework-install/include)
#include "SFibersIdentification.h"
#include "SSiPMHit.h"
#include "SSiPMCluster.h"
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


//----- scatterer geometry 4to1 (here fibers are actually SiPMs)
#define N_MODULES 1
#define N_LAYERS_PER_MODULE 4
#define N_FIBERS_PER_LAYER 28
#define N_SIDES 2

//----- canvas size
#define XCANVAS 1600
#define YCANVAS 1000

// //TOFPET 
#define NBINS_QDC 200
#define XLOW_QDC 0
#define XUP_QDC 100

#define NBINS_T 1000
#define XLOW_T 0
#define XUP_T 1000

#define NBINS_TDIFF 100
#define XLOW_TDIFF 0
#define XUP_TDIFF 10000
#define XUP_TDIFF_ZOOM 100

#define NBINS_MLR 100
#define XLOW_MLR -15
#define XUP_MLR 15

//----- cuts
#define TDIFF_MAX 10
#define QDC_MIN 0
#define T_MIN 0



/*! \file rawDisplaySiPM.C
This macro reads a sifi_results.root file (which is an output file from the sifi-framework) and produces histograms for all active addresses (module, layer, fiber) of charge and time for left and right prototype sides, as well as summed charge, average charge, time difference, MLR, correlation plots: charge left vs. charge right and time left vs. time right and a fiber multiplicity histogram. The histograms are then saved in /path/to/data/folder/YYYY_MM_DD_HH_MM/sifi_results_HISTOS.root and in /path/to/data/folder/YYYY_MM_DD_HH_MM/../pdf/2022_04_13_13_55_sifi_results_HISTOS.pdf
 */ 

void format_h_1D(TH1 * h)
{
    h->GetYaxis()->SetMaxDigits(3);
    gPad->SetLeftMargin(0.15);
    gPad->SetTopMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetBottomMargin(0.15);
}



TString ConstructPdfName(TString pathname){
    
//   TString out = pathname(0,36); // length of this path: /scratch1/gccb/data/TOFPET2/results/ - 36 characters
  TString out = pathname(0,42); // length of this path: /scratch1/gccb/data/Jan2023Beam/results_1/ - 42 characters
  out.Append("pdf/");
  TObjArray* Strings = pathname.Tokenize("/");
  TIter iString(Strings);
  TObjString* os = 0;
  TString tmp = "";
  
  while((os=(TObjString*)iString()))
  {
    tmp = os->GetString();
    if((tmp.Contains("2021") || tmp.Contains("2022") ) && tmp.Contains("_"))
      out.Append(tmp);
    else if(tmp.Contains(".root"))
    {
      out.Append("_");
      out.Append(tmp);
    }
  }
  
  out.ReplaceAll(".root",".pdf");
  delete Strings;
  
  return out;
}

struct Address {
    int iMod, iLay, iFib;
}; 

bool compareMLF(const Address &a, const Address &b)
{
    int c1, c2;
    c1 = a.iMod*N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER+a.iLay*N_FIBERS_PER_LAYER+a.iFib;
    c2 = b.iMod*N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER+b.iLay*N_FIBERS_PER_LAYER+b.iFib;
    return c1 < c2;
}
// 
Long64_t GetBasketSize(TBranch * b, bool ondisk, bool inclusive);

Long64_t GetBasketSize(TObjArray * branches, bool ondisk, bool inclusive) {
    Long64_t result = 0;
    size_t n = branches->GetEntries();
    for( size_t i = 0; i < n; ++ i ) {
       result += GetBasketSize( dynamic_cast<TBranch*>( branches->At( i ) ), ondisk, inclusive );
    }
    return result;
}

Long64_t GetBasketSize(TBranch * b, bool ondisk, bool inclusive) {
    Long64_t result = 0;
    if (b) {
       if (ondisk && b->GetZipBytes() > 0) {
          result = b->GetZipBytes();
       } else {
          result = b->GetTotBytes();
       }
       if (inclusive) {
          result += GetBasketSize(b->GetListOfBranches(), ondisk, true);
       }
       return result;
    }
    return result;
}

Long64_t GetTotalSize(TTree *t, bool ondisk) {
    TKey *key = 0;
    if (t->GetDirectory()) {
       key = t->GetDirectory()->GetKey(t->GetName());
    }
    Long64_t ondiskSize = 0;
    Long64_t totalSize = 0;
    if (key) {
       ondiskSize = key->GetNbytes();
       totalSize = key->GetObjlen();
    } else {
       TMemFile f("buffer","CREATE");
       if (t->GetCurrentFile()) {
          f.SetCompressionSettings(t->GetCurrentFile()->GetCompressionSettings());
       }
       f.WriteTObject(t);
       key = f.GetKey(t->GetName());
       ondiskSize = key->GetNbytes();
       totalSize = key->GetObjlen();
    }
    if (t->GetBranchRef() ) {
       if (ondisk) {
          ondiskSize += GetBasketSize(t->GetBranchRef(), true, true);
       } else {
          totalSize += GetBasketSize(t->GetBranchRef(), false, true);
       }
    }
    if (ondisk) {
       return ondiskSize + GetBasketSize(t->GetListOfBranches(), /* ondisk */ true, /* inclusive */ true);
    } else {
       return totalSize + GetBasketSize(t->GetListOfBranches(), /* ondisk */ false, /* inclusive */ true);
    }
}

int rawDisplaySiPMSerial(std::vector<TString> path)
{
    auto start = std::chrono::system_clock::now();
    int ii=0;
    if(!path[ii].Contains("/") || !path[ii].BeginsWith("/"))
    {
        std::cout << "##### Error! The functions needs the filename including the full absolute path..." << std::endl;
        std::abort();
    }
    gStyle->SetPalette(kBird);    
    TFile * input_file;    ///< data input file
    TTree *t;
    input_file = new TFile(path[ii]);
    if(!input_file->IsOpen()) {
        std::cerr << "##### Error in STPSource::open()! Could not open input file!" << std::endl;
        std::cerr << path[ii] << std::endl;
        return false;
    }
    t = (TTree*)input_file->Get("S");

    std::vector<std::string> str_path;
    for(int i=0; i < path.size(); i++){
        str_path.push_back(std::string(path[i]));
        printf("HERE: %s", str_path[i].c_str());
    }
    
	SLoop * loop = new SLoop();

	loop->addFiles(str_path);
	loop->setInput({});
    SCategory * pCatSiPM = SCategoryManager::getCategory(SCategory::CatSiPMHit); 
    SCategory * pCatCluster = SCategoryManager::getCategory(SCategory::CatSiPMClus); 
    vector<Address> activeAddresses;
	TH1D * hQ[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES]={nullptr};
    TH1D * hT[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES]={nullptr};

	TH1D * hSiPMMult = new TH1D("hSiPMMult", "hSiPMMult", N_FIBERS_PER_LAYER+1, -0.5, N_FIBERS_PER_LAYER+0.5); 
    TH1D * hSiPMMultPerCluster = new TH1D("hSiPMMultPerCluster", "hSiPMMultPerCluster", N_FIBERS_PER_LAYER+1, -0.5, N_FIBERS_PER_LAYER+0.5); 
    TH1D * hClusterMult = new TH1D("hClusterMult", "hClusterMult", N_FIBERS_PER_LAYER+1, -0.5, N_FIBERS_PER_LAYER+0.5); 
    
    TH2D * hSiPMHeatMapBottom = new TH2D("hSiPMHeatMapBottom", "hSiPMHeatMapBottom", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hSiPMHeatMapTop = new TH2D("hSiPMHeatMapTop", "hSiPMHeatMapTop", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    
    TH2D *hCorrelationSiPMs = new TH2D("hCorrelationSiPMs", "hCorrelationSiPMs",  N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER, 368, N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER+368,N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER, 0, N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER);
    
    TH1D * hEventTimeDiff = new TH1D("hEventTimeDiff", "hEventTimeDiff", NBINS_TDIFF, XLOW_TDIFF, XUP_TDIFF);
    TH1D * hEventTimeDiffZoom = new TH1D("hEventTimeDiffZoom", "hEventTimeDiffZoom", NBINS_TDIFF, XLOW_TDIFF, XUP_TDIFF_ZOOM);
    std::cout << "defined histograms" << std::endl;
    
    float channelQ[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES];
    
    for(int d = 0; d < N_MODULES; d++)
    {
        for(int c = 0; c < N_LAYERS_PER_MODULE; c++)
        {
            for( int a = 0; a < N_FIBERS_PER_LAYER; a++) // filling with -100
            {
                for (int b = 0; b < N_SIDES; b++)
                {
                    channelQ[d][c][a][b] = -100;
                }
            }
        }
    }
    
	Int_t mod, lay, element, cluster;
    char side;
	SSiPMHit * pHit;
    SSiPMCluster * pClus;
	UShort_t mult = 0;
    int treshold = 0;
    Long64_t time_current = 0;
    Long64_t time_prev = 0;
    Long64_t time_offset = 0;
    Long64_t ns_to_s = 10e9;
    Int_t * activeSiPMs[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES]={nullptr};

 	int nLoop = loop->getEntries();

//----- histogram definition and filling 
        loop->getEvent(0);
        pHit = (SSiPMHit *)pCatSiPM->getObject(0);
        pHit->getAddress(mod, lay, element, side);    
        time_offset = pHit->getTime();
        loop->getEvent(0);
        
	for (int i = 0; i < nLoop; ++i)
	{
		size_t nCat = pCatSiPM->getEntries();
//         std::cout << i <<" th event out of " << nLoop << std::endl;

		for (uint j = 0; j < nCat; ++j)
		{
			pHit = (SSiPMHit *)pCatSiPM->getObject(j);
			pHit->getAddress(mod, lay, element, side);
//             std::cout << mod << " " << lay << " " << element << " " << side << std::endl;
            if(mod >= N_MODULES || lay >= N_LAYERS_PER_MODULE || element >= N_FIBERS_PER_LAYER) continue;

            
            time_current = pHit->getTime();
            if(!hQ[mod][lay][element][0]){
                std::string suffix = std::string("_M") + std::to_string(mod) + std::string("L") + std::to_string(lay) + std::string("F") + std::to_string(element);
				std::string Q_L = std::string("Q") + suffix + std::string("L");
				std::string Q_R = std::string("Q") + suffix + std::string("R");
				std::string T_L = std::string("T") + suffix + std::string("L");
				std::string T_R = std::string("T") + suffix + std::string("R");

				std::string title_Q_L = Q_L + std::string("; Q_{L} [a.u.]; Counts");
				std::string title_Q_R = Q_R + std::string("; Q_{R} [a.u.]; Counts");
				std::string title_T_L = T_L + std::string("; T_{L} [s]; Counts");
				std::string title_T_R = T_R + std::string("; T_{R} [s]; Counts");
                
				hQ[mod][lay][element][0] = new TH1D(Q_L.c_str(), title_Q_L.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);
				hQ[mod][lay][element][1] = new TH1D(Q_R.c_str(), title_Q_R.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);
				hT[mod][lay][element][0] = new TH1D(T_L.c_str(), title_T_L.c_str(), NBINS_T, XLOW_T, XUP_T);
				hT[mod][lay][element][1] = new TH1D(T_R.c_str(), title_T_R.c_str(), NBINS_T, XLOW_T, XUP_T);
				
                activeAddresses.push_back({mod, lay, element});
            }     
        if(side=='l'){
            hQ[mod][lay][element][0]->Fill(pHit->getQDC());
            hT[mod][lay][element][0]->Fill((pHit->getTime()-time_offset)/ns_to_s);
            hSiPMHeatMapBottom->Fill(element, lay);
            channelQ[mod][lay][element][0]=pHit->getQDC();

        }
        if(side=='r'){      
            hQ[mod][lay][element][1]->Fill(pHit->getQDC());
            hT[mod][lay][element][1]->Fill((pHit->getTime()-time_offset)/ns_to_s);
            hSiPMHeatMapTop->Fill(element, lay);
            channelQ[mod][lay][element][1]=pHit->getQDC();
        }
                
        ++mult;
        }

          ///////////////////////////////////////////////

                SSiPMHit* pHitRight=NULL;
                SSiPMHit* pHitLeft=NULL;
                int lswID, rswID;
                //              size_t nCat = pCatSiPM->getEntries();
                for (uint right = 0; right < nCat; ++right){
                  for (uint left = 0; left < nCat; ++left){
                      pHitRight = (SSiPMHit *)pCatSiPM->getObject(right);
                      pHitRight->getAddress(mod, lay, element, side);
                      if(side!='r')  continue;

                      pHitLeft = (SSiPMHit *)pCatSiPM->getObject(left);
                      pHitLeft->getAddress(mod, lay, element, side);
                      if(side!='l')  continue;

                      int lswID, rswIF;
                      pHitLeft->getChannel(lswID);
//                       std::cout << lswID << " chr " << rswID << std::endl;
                      pHitRight->getChannel(rswID);
                      hCorrelationSiPMs->Fill(rswID,lswID);
                    }
                  }
                //      //////////////////////////////////////////////////
//     for(int kk = 0; kk < N_LAYERS_PER_MODULE; kk++)
//         {
//         for(int k = 0; k < N_FIBERS_PER_LAYER; k++) // N_LAYERS_PER_MODULE = number of fibers, this loop looks for coincidences
//         {
//             for( int ll = 0; ll<N_LAYERS_PER_MODULE; ll++) 
//             {
//             for( int l = 0; l<N_FIBERS_PER_LAYER; l++) // for every ch VS every ch
//             {
//                 if(channelQ[0][kk][k][0] > treshold && channelQ[0][ll][l][1] > treshold)
//                 {
//                     hCorrelationSiPMs->Fill(ll*N_FIBERS_PER_LAYER+l,kk*N_FIBERS_PER_LAYER+k);
//                 }
//             }
//             }
//         }
//         }
//         
//         for(int m = 0; m < N_MODULES; m++)
//         {
//             for(int l = 0; l < N_LAYERS_PER_MODULE; l++)
//             {
//                 for(int f = 0; f<N_FIBERS_PER_LAYER; f++) 
//                 {
//                     for(int s = 0; s<N_SIDES; s++) 
//                         {
//                             channelQ[m][l][f][s] = -100;
//                         }
//                 }
//             }
//         }
        hEventTimeDiff->Fill(time_current-time_prev);
        hEventTimeDiffZoom->Fill(time_current-time_prev);
        time_prev = time_current;
		hSiPMMult->Fill(mult);
		mult = 0;
		loop->nextEvent();
	}
	
	std::cout << "\n\nLoop entries: " << nLoop << std::endl;



    
    
    loop->getEvent(0);
    int mult_per_clus=0;
    int clus_mult=0;
    std::vector<Int_t> hits;
//----- histogram definition and filling    
	for (int i = 0; i < nLoop; ++i)
	{
		size_t nCatClus = pCatCluster->getEntries();
//         pCatCluster->print();
        
		for (uint j = 0; j < nCatClus; ++j)
		{
        pClus = (SSiPMCluster *)pCatCluster->getObject(j);
        if(!pClus) std::cout << "LACK OF pClus!!!!!!!" << std::endl;
//         pClus->print();
            pClus->getAddress(cluster);
            hits = pClus->getHitsArray();
            for(int k=0; k< hits.size(); k++)
            {
                mult_per_clus++;
//                 hSiPMMultPerCluster->Fill(mult_per_clus); //if hSiPMMultPerCluster is filled here, we obtain the contents of the SSiPMCluster.data.hits branch
            }
            hSiPMMultPerCluster->Fill(mult_per_clus);
            mult_per_clus=0;
            clus_mult++;
//             hClusterMult->Fill(clus_mult); //if hClusterMult filled here, we obtain the contents of the SSiPMCluster.data.cluster branch
        }
        hClusterMult->Fill(clus_mult);
        clus_mult=0;
        loop->nextEvent();
    }
    
    
	path[ii].ReplaceAll(".root", "_HISTOS.root");
	TFile *output = new TFile(path[ii],"RECREATE");
	output->cd();

//----- canvas definition    
    TCanvas * canQDC_layer0 = new TCanvas("canQDC_layer0", "canQDC_layer0", XCANVAS, YCANVAS); 
    canQDC_layer0->SetLeftMargin(0.2);
	canQDC_layer0->DivideSquare(activeAddresses.size()/4+1);

    TCanvas * canQDC_layer1 = new TCanvas("canQDC_layer1", "canQDC_layer1", XCANVAS, YCANVAS); 
    canQDC_layer1->SetLeftMargin(0.2);
	canQDC_layer1->DivideSquare(activeAddresses.size()/4+1);
    
    TCanvas * canQDC_layer2 = new TCanvas("canQDC_layer2", "canQDC_layer2", XCANVAS, YCANVAS); 
    canQDC_layer2->SetLeftMargin(0.2);
	canQDC_layer2->DivideSquare(activeAddresses.size()/4+1);
    
    TCanvas * canQDC_layer3 = new TCanvas("canQDC_layer3", "canQDC_layer3", XCANVAS, YCANVAS); 
    canQDC_layer3->SetLeftMargin(0.2);
	canQDC_layer3->DivideSquare(activeAddresses.size()/4+1);
    
    
    TCanvas * canTime_layer0 = new TCanvas("canTime_layer0", "canTime_layer0", XCANVAS, YCANVAS); 
    canTime_layer0->SetLeftMargin(0.2);
	canTime_layer0->DivideSquare(activeAddresses.size()/4+1);

    TCanvas * canTime_layer1 = new TCanvas("canTime_layer1", "canTime_layer1", XCANVAS, YCANVAS); 
    canTime_layer1->SetLeftMargin(0.2);
	canTime_layer1->DivideSquare(activeAddresses.size()/4+1);
    
    TCanvas * canTime_layer2 = new TCanvas("canTime_layer2", "canTime_layer2", XCANVAS, YCANVAS); 
    canTime_layer2->SetLeftMargin(0.2);
	canTime_layer2->DivideSquare(activeAddresses.size()/4+1);
    
    TCanvas * canTime_layer3 = new TCanvas("canTime_layer3", "canTime_layer3", XCANVAS, YCANVAS); 
    canTime_layer3->SetLeftMargin(0.2);
	canTime_layer3->DivideSquare(activeAddresses.size()/4+1);
    
    
	TCanvas * canMult = new TCanvas("Mult", "Mult", 600, 400); 
    canMult->SetLeftMargin(0.2);
    
    TCanvas * canSiPMMultPerCluster = new TCanvas("SiPMMultPerCluster", "SiPMMultPerCluster", 600, 400); 
    canSiPMMultPerCluster->SetLeftMargin(0.2);
    
    TCanvas * canClusterMult = new TCanvas("ClusterMult", "ClusterMult", 600, 400); 
    canClusterMult->SetLeftMargin(0.2);
    
    TCanvas * canSiPMHeatMap = new TCanvas("SiPMHeatMap", "SiPMHeatMap", 1200, 400);
    canSiPMHeatMap->SetLeftMargin(0.2);
    canSiPMHeatMap->Divide(2,1);
    
    TCanvas * canSiPMCorrelations = new TCanvas("SiPMCorrelations", "SiPMCorrelations", 600, 400);
    canSiPMCorrelations->SetLeftMargin(0.2);
    
    TCanvas * canEventTimeDiff = new TCanvas("EventTimeDiff", "EventTimeDiff", 1200, 400);
    canEventTimeDiff->SetLeftMargin(0.2);
    canEventTimeDiff->Divide(2,1);
    
    Int_t m=0;
    Int_t l=0;
    Int_t f=0;
    Int_t colL = kBlack;
    Int_t colR = kRed;

std::sort(activeAddresses.begin(), activeAddresses.end(), compareMLF);    

//----- canvas plotting  

    for(int i=0; i< activeAddresses.size()/4; i++){
        m=activeAddresses[i].iMod;
        l=activeAddresses[i].iLay;
        f=activeAddresses[i].iFib;

        if(!output->GetListOfKeys()->Contains(Form("Module%d", m))) output->mkdir(Form("Module%d/", m));
        output->cd(Form("Module%d/", m));
        if(!gDirectory->GetListOfKeys()->Contains(Form("Layer%d", l))) output->mkdir(Form("Module%d/Layer%d", m, l));
        output->cd(Form("Module%d/Layer%d", m, l));
        if(!gDirectory->GetListOfKeys()->Contains(Form("Fiber%d", f))) output->mkdir(Form("Module%d/Layer%d/Fiber%d",m, l, f));
        output->cd(Form("Module%d/Layer%d/Fiber%d",m, l, f));  
        hQ[m][l][f][0]->SetDirectory(gDirectory);
        hQ[m][l][f][0]->Write();
        hQ[m][l][f][1]->SetDirectory(gDirectory);
        hQ[m][l][f][1]->Write();
        hT[m][l][f][0]->SetDirectory(gDirectory);
        hT[m][l][f][0]->Write();
        hT[m][l][f][1]->SetDirectory(gDirectory);
        hT[m][l][f][1]->Write();
        
        canQDC_layer0->cd(i+1);
        gPad->SetGrid(1, 1);
        hQ[m][l][f][0]->SetLineColor(colL);
        hQ[m][l][f][1]->SetLineColor(colR);
        hQ[m][l][f][0]->Draw();
        hQ[m][l][f][1]->Draw("same");
        hQ[m][l][f][0]->GetYaxis()->SetRangeUser(0.1, TMath::Max(hQ[m][l][f][0]->GetBinContent(hQ[m][l][f][0]->GetMaximumBin()), hQ[m][l][f][1]->GetBinContent(hQ[m][l][f][1]->GetMaximumBin())) + 50);
        gPad->SetLogy();
        format_h_1D(hQ[m][l][f][0]);
        format_h_1D(hQ[m][l][f][1]);
        
        canTime_layer0->cd(i+1);
        gPad->SetGrid(1, 1);
        hT[m][l][f][0]->SetLineColor(colL);
        hT[m][l][f][1]->SetLineColor(colR);
        hT[m][l][f][0]->GetXaxis()->SetRangeUser(hT[m][l][f][0]->GetXaxis()->GetBinCenter(hT[m][l][f][0]->FindFirstBinAbove())-1, hT[m][l][f][0]->GetXaxis()->GetBinCenter(hT[m][l][f][0]->FindLastBinAbove())+1);
        hT[m][l][f][0]->Draw();
        hT[m][l][f][1]->Draw("same");
        hT[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hT[m][l][f][0]->GetBinContent(hT[m][l][f][0]->GetMaximumBin()), hT[m][l][f][1]->GetBinContent(hT[m][l][f][1]->GetMaximumBin())) + 50);
        format_h_1D(hT[m][l][f][0]);
        format_h_1D(hT[m][l][f][1]); 
        }
        
    for(int i = activeAddresses.size()/4; i < activeAddresses.size()/2; i++){
        m=activeAddresses[i].iMod;
        l=activeAddresses[i].iLay;
        f=activeAddresses[i].iFib;
        
        if(!output->GetListOfKeys()->Contains(Form("Module%d", m))) output->mkdir(Form("Module%d/", m));
        output->cd(Form("Module%d/", m));
        if(!gDirectory->GetListOfKeys()->Contains(Form("Layer%d", l))) output->mkdir(Form("Module%d/Layer%d", m, l));
        output->cd(Form("Module%d/Layer%d", m, l));
        if(!gDirectory->GetListOfKeys()->Contains(Form("Fiber%d", f))) output->mkdir(Form("Module%d/Layer%d/Fiber%d",m, l, f));
        output->cd(Form("Module%d/Layer%d/Fiber%d",m, l, f));  
        hQ[m][l][f][0]->SetDirectory(gDirectory);
        hQ[m][l][f][0]->Write();
        hQ[m][l][f][1]->SetDirectory(gDirectory);
        hQ[m][l][f][1]->Write();
        hT[m][l][f][0]->SetDirectory(gDirectory);
        hT[m][l][f][0]->Write();
        hT[m][l][f][1]->SetDirectory(gDirectory);
        hT[m][l][f][1]->Write();
        
        canQDC_layer1->cd(i-activeAddresses.size()/4+1);
        gPad->SetGrid(1, 1);
        hQ[m][l][f][0]->SetLineColor(colL);
        hQ[m][l][f][1]->SetLineColor(colR);
        hQ[m][l][f][0]->Draw();
        hQ[m][l][f][1]->Draw("same");
        hQ[m][l][f][0]->GetYaxis()->SetRangeUser(0.1, TMath::Max(hQ[m][l][f][0]->GetBinContent(hQ[m][l][f][0]->GetMaximumBin()), hQ[m][l][f][1]->GetBinContent(hQ[m][l][f][1]->GetMaximumBin())) + 50);
        gPad->SetLogy();
        format_h_1D(hQ[m][l][f][0]);
        format_h_1D(hQ[m][l][f][1]);
        
        canTime_layer1->cd(i-activeAddresses.size()/4+1);
        gPad->SetGrid(1, 1);
        hT[m][l][f][0]->SetLineColor(colL);
        hT[m][l][f][1]->SetLineColor(colR);
        hT[m][l][f][0]->GetXaxis()->SetRangeUser(hT[m][l][f][0]->GetXaxis()->GetBinCenter(hT[m][l][f][0]->FindFirstBinAbove())-1, hT[m][l][f][0]->GetXaxis()->GetBinCenter(hT[m][l][f][0]->FindLastBinAbove())+1);
        hT[m][l][f][0]->Draw();
        hT[m][l][f][1]->Draw("same");
        hT[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hT[m][l][f][0]->GetBinContent(hT[m][l][f][0]->GetMaximumBin()), hT[m][l][f][1]->GetBinContent(hT[m][l][f][1]->GetMaximumBin())) + 50);
        format_h_1D(hT[m][l][f][0]);
        format_h_1D(hT[m][l][f][1]);
        }
        
    for(int i = activeAddresses.size()/2; i < 3*activeAddresses.size()/4; i++){
        m=activeAddresses[i].iMod;
        l=activeAddresses[i].iLay;
        f=activeAddresses[i].iFib;
        
        if(!output->GetListOfKeys()->Contains(Form("Module%d", m))) output->mkdir(Form("Module%d/", m));
        output->cd(Form("Module%d/", m));
        if(!gDirectory->GetListOfKeys()->Contains(Form("Layer%d", l))) output->mkdir(Form("Module%d/Layer%d", m, l));
        output->cd(Form("Module%d/Layer%d", m, l));
        if(!gDirectory->GetListOfKeys()->Contains(Form("Fiber%d", f))) output->mkdir(Form("Module%d/Layer%d/Fiber%d",m, l, f));
        output->cd(Form("Module%d/Layer%d/Fiber%d",m, l, f));  
        hQ[m][l][f][0]->SetDirectory(gDirectory);
        hQ[m][l][f][0]->Write();
        hQ[m][l][f][1]->SetDirectory(gDirectory);
        hQ[m][l][f][1]->Write();
        hT[m][l][f][0]->SetDirectory(gDirectory);
        hT[m][l][f][0]->Write();
        hT[m][l][f][1]->SetDirectory(gDirectory);
        hT[m][l][f][1]->Write();
        
        canQDC_layer2->cd(i-activeAddresses.size()/2+1);
        gPad->SetGrid(1, 1);
        hQ[m][l][f][0]->SetLineColor(colL);
        hQ[m][l][f][1]->SetLineColor(colR);
        hQ[m][l][f][0]->Draw();
        hQ[m][l][f][1]->Draw("same");
        hQ[m][l][f][0]->GetYaxis()->SetRangeUser(0.1, TMath::Max(hQ[m][l][f][0]->GetBinContent(hQ[m][l][f][0]->GetMaximumBin()), hQ[m][l][f][1]->GetBinContent(hQ[m][l][f][1]->GetMaximumBin())) + 50);
        gPad->SetLogy();
        format_h_1D(hQ[m][l][f][0]);
        format_h_1D(hQ[m][l][f][1]);
        
        canTime_layer2->cd(i-activeAddresses.size()/2+1);
        gPad->SetGrid(1, 1);
        hT[m][l][f][0]->SetLineColor(colL);
        hT[m][l][f][1]->SetLineColor(colR);
        hT[m][l][f][0]->GetXaxis()->SetRangeUser(hT[m][l][f][0]->GetXaxis()->GetBinCenter(hT[m][l][f][0]->FindFirstBinAbove())-1, hT[m][l][f][0]->GetXaxis()->GetBinCenter(hT[m][l][f][0]->FindLastBinAbove())+1);
        hT[m][l][f][0]->Draw();
        hT[m][l][f][1]->Draw("same");
        hT[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hT[m][l][f][0]->GetBinContent(hT[m][l][f][0]->GetMaximumBin()), hT[m][l][f][1]->GetBinContent(hT[m][l][f][1]->GetMaximumBin())) + 50);
        format_h_1D(hT[m][l][f][0]);
        format_h_1D(hT[m][l][f][1]);
        }
        
    for(int i = 3*activeAddresses.size()/4; i < activeAddresses.size(); i++){
        m=activeAddresses[i].iMod;
        l=activeAddresses[i].iLay;
        f=activeAddresses[i].iFib;
        
        if(!output->GetListOfKeys()->Contains(Form("Module%d", m))) output->mkdir(Form("Module%d/", m));
        output->cd(Form("Module%d/", m));
        if(!gDirectory->GetListOfKeys()->Contains(Form("Layer%d", l))) output->mkdir(Form("Module%d/Layer%d", m, l));
        output->cd(Form("Module%d/Layer%d", m, l));
        if(!gDirectory->GetListOfKeys()->Contains(Form("Fiber%d", f))) output->mkdir(Form("Module%d/Layer%d/Fiber%d",m, l, f));
        output->cd(Form("Module%d/Layer%d/Fiber%d",m, l, f));  
        hQ[m][l][f][0]->SetDirectory(gDirectory);
        hQ[m][l][f][0]->Write();
        hQ[m][l][f][1]->SetDirectory(gDirectory);
        hQ[m][l][f][1]->Write();
        hT[m][l][f][0]->SetDirectory(gDirectory);
        hT[m][l][f][0]->Write();
        hT[m][l][f][1]->SetDirectory(gDirectory);
        hT[m][l][f][1]->Write();
        
        canQDC_layer3->cd(i-3*activeAddresses.size()/4+1);
        gPad->SetGrid(1, 1);
        hQ[m][l][f][0]->SetLineColor(colL);
        hQ[m][l][f][1]->SetLineColor(colR);
        hQ[m][l][f][0]->Draw();
        hQ[m][l][f][1]->Draw("same");
        hQ[m][l][f][0]->GetYaxis()->SetRangeUser(0.1, TMath::Max(hQ[m][l][f][0]->GetBinContent(hQ[m][l][f][0]->GetMaximumBin()), hQ[m][l][f][1]->GetBinContent(hQ[m][l][f][1]->GetMaximumBin())) + 50);
         gPad->SetLogy();
        format_h_1D(hQ[m][l][f][0]);
        format_h_1D(hQ[m][l][f][1]);
        
        canTime_layer3->cd(i-3*activeAddresses.size()/4+1);
        gPad->SetGrid(1, 1);
        hT[m][l][f][0]->SetLineColor(colL);
        hT[m][l][f][1]->SetLineColor(colR);
        hT[m][l][f][0]->GetXaxis()->SetRangeUser(hT[m][l][f][0]->GetXaxis()->GetBinCenter(hT[m][l][f][0]->FindFirstBinAbove())-1, hT[m][l][f][0]->GetXaxis()->GetBinCenter(hT[m][l][f][0]->FindLastBinAbove())+1);
        hT[m][l][f][0]->Draw();
        hT[m][l][f][1]->Draw("same");
        hT[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hT[m][l][f][0]->GetBinContent(hT[m][l][f][0]->GetMaximumBin()), hT[m][l][f][1]->GetBinContent(hT[m][l][f][1]->GetMaximumBin())) + 50);
        format_h_1D(hT[m][l][f][0]);
        format_h_1D(hT[m][l][f][1]);
        }

	canMult->cd();  
    gPad->SetGrid(1, 1);
	hSiPMMult->Draw();
 	format_h_1D(hSiPMMult);
    hSiPMMult->GetXaxis()->SetRange(1, activeAddresses.size()+1);
	hSiPMMult->GetXaxis()->SetTitle("# SiPMs");
	hSiPMMult->GetYaxis()->SetTitle("Counts");
    
    canSiPMMultPerCluster->cd();  
    gPad->SetGrid(1, 1);
	hSiPMMultPerCluster->Draw();
 	format_h_1D(hSiPMMultPerCluster);
    hSiPMMultPerCluster->GetXaxis()->SetRange(1, activeAddresses.size()+1);
	hSiPMMultPerCluster->GetXaxis()->SetTitle("# SiPMs per cluster");
	hSiPMMultPerCluster->GetYaxis()->SetTitle("Counts");

    canClusterMult->cd();  
    gPad->SetGrid(1, 1);
	hClusterMult->Draw();
 	format_h_1D(hClusterMult);
    hClusterMult->GetXaxis()->SetRange(1, activeAddresses.size()+1);
    hClusterMult->GetXaxis()->SetTitle("# clusters per event");
	hClusterMult->GetYaxis()->SetTitle("Counts");
    
    
    
    canSiPMHeatMap->cd(1);
    gPad->SetGrid(1, 1);
    hSiPMHeatMapBottom->Draw("colz");
    hSiPMHeatMapBottom->GetXaxis()->SetTitle("SiPM # in a layer (\"element\")");
    hSiPMHeatMapBottom->GetYaxis()->SetTitle("Layer #");
      
    canSiPMHeatMap->cd(2);
    gPad->SetGrid(1, 1);
    hSiPMHeatMapTop->Draw("colz");
    hSiPMHeatMapTop->GetXaxis()->SetTitle("SiPM # in a layer (\"element\")");
    hSiPMHeatMapTop->GetYaxis()->SetTitle("Layer #");
    
    canSiPMCorrelations->cd();
    hCorrelationSiPMs->GetXaxis()->SetTitle("Right ch");
    hCorrelationSiPMs->GetYaxis()->SetTitle("Left ch");
    hCorrelationSiPMs->GetXaxis()->SetNdivisions(004);
    hCorrelationSiPMs->GetYaxis()->SetNdivisions(004);
    hCorrelationSiPMs->Draw("COLZ");
    TLine *l1 = new TLine(0,0,N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER,N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER);
    TLine *l2 = new TLine(0,N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER,N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER,0);
    l1->Draw();
    l2->Draw();
    
    canEventTimeDiff->cd();
    canEventTimeDiff->cd(1);
    gPad->SetGrid(1, 1);
    gPad->SetLogy();
    hEventTimeDiff->Draw();
    hEventTimeDiff->GetXaxis()->SetTitle("t_current - t_prev [ns]");
    hEventTimeDiff->GetYaxis()->SetTitle("Counts");
    
    canEventTimeDiff->cd(2);
    gPad->SetGrid(1, 1);
    gPad->SetLogy();
    hEventTimeDiffZoom->Draw();
    hEventTimeDiffZoom->GetXaxis()->SetTitle("t_current - t_prev [ns]");
    hEventTimeDiffZoom->GetYaxis()->SetTitle("Counts");
    
	output->mkdir("Canvases");
	output->cd("Canvases");
    
	canQDC_layer0->Write();
    canQDC_layer1->Write();
    canQDC_layer2->Write();
    canQDC_layer3->Write();
    canTime_layer0->Write();
    canTime_layer1->Write();
    canTime_layer2->Write();
    canTime_layer3->Write();
    canMult->Write();
    canSiPMMultPerCluster->Write();
    canClusterMult->Write();
    canSiPMHeatMap->Write();
    canSiPMCorrelations->Write();
    canEventTimeDiff->Write();
//     output->Close();

// //----- writing canvas to file    
    TString pdfname = ConstructPdfName(path[ii]);
    canQDC_layer0->Print(Form("%s(",pdfname.Data()), "pdf");
    canQDC_layer1->Print(pdfname,"pdf");
    canQDC_layer2->Print(pdfname,"pdf");
    canQDC_layer3->Print(pdfname,"pdf");
    canTime_layer0->Print(pdfname,"pdf");
    canTime_layer1->Print(pdfname,"pdf");
    canTime_layer2->Print(pdfname,"pdf");
    canTime_layer3->Print(pdfname,"pdf");
    canMult->Print(pdfname,"pdf");
    canSiPMMultPerCluster->Print(pdfname,"pdf");
    canClusterMult->Print(pdfname,"pdf");
    canSiPMHeatMap->Print(pdfname,"pdf");
    canSiPMCorrelations->Print(pdfname,"pdf");
    canEventTimeDiff->Print(Form("%s)",pdfname.Data()),"pdf");
    
    auto end = std::chrono::system_clock::now();   
    std::chrono::duration<double> time = end-start;
    std::cout << "Time: " << time.count() << " s\n";
    std::cout << "Tree size: " << GetTotalSize(t, true)/1000000. << " MB\n"; 
    std::cout << "Entries: " << nLoop << "\n";
    std::cout << "Processing time per loop entry: " << time.count()/nLoop << " s/entry"<< std::endl;
    std::cout << "Processing time per MB: " << time.count()/GetTotalSize(t, true)*1000000. << " s/MB"<< std::endl;
    
	return 0;
}
