#include "SCategoryManager.h"
#include "SLoop.h"
#include "SDDSamples.h"
#include "SFibersRaw.h"
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


//----- prototype geometry 4to1
#define N_MODULES 1
#define N_LAYERS_PER_MODULE 7
#define N_FIBERS_PER_LAYER 55
#define N_SIDES 2

//----- canvas size
#define XCANVAS 1600
#define YCANVAS 1000

// //TOFPET 
#define NBINS_QDC 200
#define XLOW_QDC 0
#define XUP_QDC 100

#define NBINS_T 100
#define XLOW_T 0
#define XUP_T 5e10

#define NBINS_TDIFF 100
#define XLOW_TDIFF -50
#define XUP_TDIFF 50

#define NBINS_MLR 100
#define XLOW_MLR -15
#define XUP_MLR 15

//----- cuts
#define TDIFF_MAX 10
#define QDC_MIN 0
#define T_MIN 0


/*! \file rawDisplayDynamicDAQ.C
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




int rawDisplayDynamicDAQ(TString path)
{
    auto start = std::chrono::system_clock::now();  
    
    if(!path.Contains("/") || !path.BeginsWith("/"))
    {
        std::cout << "##### Error! The functions needs the filename including the full absolute path..." << std::endl;
        std::abort();
    }
    
    gStyle->SetPalette(kBird);
	SLoop * loop = new SLoop();
	loop->addFile(std::string(path));
	loop->setInput({});

	SCategory * pCatRaw = SCategoryManager::getCategory(SCategory::CatFibersRaw);
// 	SCategory * pCatSample = SCategoryManager::getCategory(SCategory::CatDDSamples);
    

    
    vector<Address> activeAddresses;
//     SiFiDatabase * db = new SiFiDatabase("/path/to/database.db");
//     SiFiSeries * sifi = db->GetSeries(SERIES_NO)
//     int N_MODULES = sifi->fModules;
//     for(int i = 0; i< N_MODULES; i++){
// //         sifi->fLayers;    
//     }
    
    TH1D * hA[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES]={nullptr};
	TH1D * hQ[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES]={nullptr};
    TH1D * hT[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES]={nullptr};
    
	TH1D * hQSum[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER]={nullptr};
    TH1D * hQAve[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER]={nullptr};
	TH1D * hTDif[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER]={nullptr};
    TH1D * hMLR[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER]={nullptr};
    
	TH2D * hQLvsQR[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER]={nullptr};
    TH2D * hTLvsTR[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER]={nullptr};
    
	TH1D * hFiberMult = new TH1D("hFiberMult", "hSiPMMult", N_FIBERS_PER_LAYER+1, -0.5, N_FIBERS_PER_LAYER+0.5); 

	Int_t mod, lay, fi, side;
	SFibersRaw * pRaw;
	UShort_t mult = 0;
    
    Int_t * activeSiPMs[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES]={nullptr};
    
    
    

    
    
 	int nLoop = loop->getEntries();

//----- histogram definition and filling    
	for (int i = 0; i < nLoop; ++i)
	{
		size_t nCat = pCatRaw->getEntries();
        
		for (uint j = 0; j < nCat; ++j)
		{
			pRaw = (SFibersRaw *)pCatRaw->getObject(j);
			pRaw->getAddress(mod, lay, fi);
            
            if(!hQ[mod][lay][fi][0]){
//             if(!gDirectory->FindObject(hQ[mod][lay][fi][0]))
                std::string suffix = std::string("_M") + std::to_string(mod) + std::string("L") + std::to_string(lay) + std::string("F") + std::to_string(fi);
				std::string Q_L = std::string("Q") + suffix + std::string("L");
				std::string Q_R = std::string("Q") + suffix + std::string("R");
				std::string T_L = std::string("T") + suffix + std::string("L");
				std::string T_R = std::string("T") + suffix + std::string("R");
				std::string TDif = std::string("TDif") + suffix;
                std::string QAve = std::string("QAve") + suffix;
				std::string QSum = std::string("QSum") + suffix;
				std::string MLR = std::string("MLR") + suffix;
				std::string QLvsQR = std::string("QLvsQR") + suffix;
                std::string TLvsTR = std::string("TLvsTR") + suffix;

				std::string title_Q_L = Q_L + std::string("; Q_{L} [a.u.]; Counts");
				std::string title_Q_R = Q_R + std::string("; Q_{R} [a.u.]; Counts");
				std::string title_T_L = T_L + std::string("; T_{L} [ns]; Counts");
				std::string title_T_R = T_R + std::string("; T_{R} [ns]; Counts");
				std::string title_TDif = TDif + std::string("; T_{dif} [ns]; Counts");
                std::string title_QAve = QAve + std::string("; Q_{ave} [a.u.]; Counts");
				std::string title_QSum = QSum + std::string("; Q_{L} + Q_{R} [a.u.]; Counts");
                std::string title_MLR = MLR + std::string("; M_{LR}; Counts");
				std::string title_QLvsQR = QLvsQR + std::string("; Q_{L} [a.u.]; Q_{R} [a.u.]");
				std::string title_TLvsTR = TLvsTR + std::string("; T_{L} [ns]; T_{R} [ns]");
                
				hQ[mod][lay][fi][0] = new TH1D(Q_L.c_str(), title_Q_L.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);
				hQ[mod][lay][fi][1] = new TH1D(Q_R.c_str(), title_Q_R.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);
				hT[mod][lay][fi][0] = new TH1D(T_L.c_str(), title_T_L.c_str(), NBINS_T, XLOW_T, XUP_T);
				hT[mod][lay][fi][1] = new TH1D(T_R.c_str(), title_T_R.c_str(), NBINS_T, XLOW_T, XUP_T);
				
                hTDif[mod][lay][fi] = new TH1D(TDif.c_str(), title_TDif.c_str(), NBINS_TDIFF, XLOW_TDIFF, XUP_TDIFF);
                hQAve[mod][lay][fi] = new TH1D(QAve.c_str(), title_QAve.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);
				hQSum[mod][lay][fi] = new TH1D(QSum.c_str(), title_QSum.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC*2);
                hMLR[mod][lay][fi] = new TH1D(MLR.c_str(), title_MLR.c_str(), NBINS_MLR, XLOW_MLR, XUP_MLR);
				
                hQLvsQR[mod][lay][fi] = new TH2D(QLvsQR.c_str(), title_QLvsQR.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC, NBINS_QDC, XLOW_QDC, XUP_QDC);
                hTLvsTR[mod][lay][fi] = new TH2D(TLvsTR.c_str(), title_TLvsTR.c_str(), NBINS_T, XLOW_T, XUP_T, NBINS_T, XLOW_T, XUP_T );
                
                activeAddresses.push_back({mod, lay, fi});
            }     
//             cout.precision(3);
//             std::cout << "EEE "<<j << " "  << i << "   "  << mod << " " << lay << " " << fi << "    " << pRaw->getQDCL() << " " << pRaw->getTimeL() << "    " << pRaw->getQDCR() << " " << pRaw->getTimeR() << " " << /*pRaw->getAddress(mod, lay, fi) <<*/ std::endl;
//             std::cout << "NOW" << std::endl;
//             std::cout << "left: " << pRaw->getQDCL() << std::endl;     
//             std::cout << "right: " << pRaw->getQDCR() << std::endl;  
//             std::cout << std::endl;
//             std::cout << mod << " " << lay << " " << fi << std::endl;
//             if(hQ[mod][lay][fi][0])  std::cout << "EXISTS" << std::endl;
//             hQ[mod][lay][fi][0]->Fill(pRaw->getQDCL());
//             if(pRaw->getQDCL() > QDC_MIN &&
//             pRaw->getQDCR() > QDC_MIN &&
//             pRaw->getTimeL() > T_MIN &&
//             pRaw->getTimeR() > T_MIN /*&&
//             fabs(pRaw->getTimeL() - pRaw->getTimeR()) < TDIFF_MAX*/)
//             {
                hQ[mod][lay][fi][0]->Fill(pRaw->getQDCL());
                hQ[mod][lay][fi][1]->Fill(pRaw->getQDCR());
                hT[mod][lay][fi][0]->Fill(pRaw->getTimeL());
                hT[mod][lay][fi][1]->Fill(pRaw->getTimeR());
                
                hQAve[mod][lay][fi]->Fill(sqrt(pRaw->getQDCL()*pRaw->getQDCR()));
                hQSum[mod][lay][fi]->Fill(pRaw->getQDCL() + pRaw->getQDCR());
                hMLR[mod][lay][fi]->Fill(log(sqrt(pRaw->getQDCR()/pRaw->getQDCL())));
                
                hQLvsQR[mod][lay][fi]->Fill(pRaw->getQDCL(), pRaw->getQDCR());
                hTLvsTR[mod][lay][fi]->Fill(pRaw->getTimeL(), pRaw->getTimeR());
                
//                 if(fabs(pRaw->getTimeL() - pRaw->getTimeR() ) < TDIFF_MAX)
//                 {
                    ++mult;
//                 }
                
//             }	

//             if(pRaw->getQDCL() > QDC_MIN &&
//             pRaw->getQDCR() > QDC_MIN &&
//             pRaw->getTimeL() > T_MIN &&
//             pRaw->getTimeR() > T_MIN)
//             {
                hTDif[mod][lay][fi]->Fill(pRaw->getTimeL() - pRaw->getTimeR() );
//             }
			
		}
		
		hFiberMult->Fill(mult);
		mult = 0;
		loop->nextEvent();
	}
	
	std::cout << "\n\nLoop entries: " << nLoop << std::endl;

	path.ReplaceAll(".root", "_HISTOS.root");
	TFile *output = new TFile(path,"RECREATE");
	output->cd();

//----- canvas definition    
    TCanvas * canQDC_layer0 = new TCanvas("canQDC_layer0", "canQDC_layer0", XCANVAS, YCANVAS); 
    canQDC_layer0->SetLeftMargin(0.2);
	canQDC_layer0->DivideSquare(activeAddresses.size()/4);

    TCanvas * canQDC_layer1 = new TCanvas("canQDC_layer1", "canQDC_layer1", XCANVAS, YCANVAS); 
    canQDC_layer1->SetLeftMargin(0.2);
	canQDC_layer1->DivideSquare(activeAddresses.size()/4);
    
    TCanvas * canQDC_layer2 = new TCanvas("canQDC_layer2", "canQDC_layer2", XCANVAS, YCANVAS); 
    canQDC_layer2->SetLeftMargin(0.2);
	canQDC_layer2->DivideSquare(activeAddresses.size()/4);
    
    TCanvas * canQDC_layer3 = new TCanvas("canQDC_layer3", "canQDC_layer3", XCANVAS, YCANVAS); 
    canQDC_layer3->SetLeftMargin(0.2);
	canQDC_layer3->DivideSquare(activeAddresses.size()/4);
//     TCanvas * canT = new TCanvas("T", "T", XCANVAS, YCANVAS); 
//     canT->SetLeftMargin(0.2); 
//     canT->DivideSquare(activeAddresses.size());
//     
// 	TCanvas * canTDif = new TCanvas("TDif", "TDif", XCANVAS, YCANVAS);
//     canTDif->SetLeftMargin(0.2);
// 	canTDif->DivideSquare(activeAddresses.size());
//     
//     TCanvas * canQAve = new TCanvas("QAve", "QAve", XCANVAS, YCANVAS); 
//     canQAve->SetLeftMargin(0.2); 
//     canQAve->DivideSquare(activeAddresses.size());
//     
// 	TCanvas * canQSum = new TCanvas("QSum", "QSum", XCANVAS, YCANVAS); 
//     canQSum->SetLeftMargin(0.2); 
// 	canQSum->DivideSquare(activeAddresses.size());
//     
// 	TCanvas * canMLR = new TCanvas("MLR", "MLR", XCANVAS, YCANVAS); 
//     canMLR->SetLeftMargin(0.2); 
// 	canMLR->DivideSquare(activeAddresses.size());
//     
// 	TCanvas * canQLvsQR = new TCanvas("QLvsQR", "Q_{L} vs Q_{R}", XCANVAS, YCANVAS); 
//     canQLvsQR->SetLeftMargin(0.2);
// 	canQLvsQR->DivideSquare(activeAddresses.size());
//     
//     TCanvas * canTLvsTR = new TCanvas("TLvsTR", "T_{L} vs T_{R}", XCANVAS, YCANVAS); 
//     canTLvsTR->SetLeftMargin(0.2);
// 	canTLvsTR->DivideSquare(activeAddresses.size());
    
	TCanvas * canMult = new TCanvas("Mult", "Mult", 600, 400); 
    canMult->SetLeftMargin(0.2);
    
    
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
 
//         hQ[m][l][f][0]->SetDirectory(gDirectory);
//         hQ[m][l][f][0]->Write();
// 

//         hT[m][l][f][0]->SetDirectory(gDirectory);
//         hT[m][l][f][0]->Write();
//         hT[m][l][f][1]->SetDirectory(gDirectory);
//         hT[m][l][f][1]->Write();
// 
//         hTDif[m][l][f]->SetDirectory(gDirectory);
//         hTDif[m][l][f]->Write();
//         hQSum[m][l][f]->SetDirectory(gDirectory);
//         hQSum[m][l][f]->Write();
//         hQAve[m][l][f]->SetDirectory(gDirectory);
//         hQAve[m][l][f]->Write();
//         hMLR[m][l][f]->SetDirectory(gDirectory);
//         hMLR[m][l][f]->Write();
//         
//         hQLvsQR[m][l][f]->SetDirectory(gDirectory);
//         hQLvsQR[m][l][f]->Write();
//         hTLvsTR[m][l][f]->SetDirectory(gDirectory);
//         hTLvsTR[m][l][f]->Write();
        

        canQDC_layer0->cd(i+1);
        gPad->SetGrid(1, 1);
        hQ[m][l][f][0]->SetLineColor(colL);
        hQ[m][l][f][1]->SetLineColor(colR);
        hQ[m][l][f][0]->Draw();
        hQ[m][l][f][1]->Draw("same");
        hQ[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hQ[m][l][f][0]->GetBinContent(hQ[m][l][f][0]->GetMaximumBin()), hQ[m][l][f][1]->GetBinContent(hQ[m][l][f][1]->GetMaximumBin())) + 50);
        format_h_1D(hQ[m][l][f][0]);
        format_h_1D(hQ[m][l][f][1]);
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
 
        
        canQDC_layer1->cd(i-activeAddresses.size()/4+1);
        gPad->SetGrid(1, 1);
        hQ[m][l][f][0]->SetLineColor(colL);
        hQ[m][l][f][1]->SetLineColor(colR);
        hQ[m][l][f][0]->Draw();
        hQ[m][l][f][1]->Draw("same");
        hQ[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hQ[m][l][f][0]->GetBinContent(hQ[m][l][f][0]->GetMaximumBin()), hQ[m][l][f][1]->GetBinContent(hQ[m][l][f][1]->GetMaximumBin())) + 50);
        format_h_1D(hQ[m][l][f][0]);
        format_h_1D(hQ[m][l][f][1]);
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
        
        canQDC_layer2->cd(i-activeAddresses.size()/2+1);
        gPad->SetGrid(1, 1);
        hQ[m][l][f][0]->SetLineColor(colL);
        hQ[m][l][f][1]->SetLineColor(colR);
        hQ[m][l][f][0]->Draw();
        hQ[m][l][f][1]->Draw("same");
        hQ[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hQ[m][l][f][0]->GetBinContent(hQ[m][l][f][0]->GetMaximumBin()), hQ[m][l][f][1]->GetBinContent(hQ[m][l][f][1]->GetMaximumBin())) + 50);
        format_h_1D(hQ[m][l][f][0]);
        format_h_1D(hQ[m][l][f][1]);
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
        
        canQDC_layer3->cd(i-3*activeAddresses.size()/4+1);
        gPad->SetGrid(1, 1);
        hQ[m][l][f][0]->SetLineColor(colL);
        hQ[m][l][f][1]->SetLineColor(colR);
        hQ[m][l][f][0]->Draw();
        hQ[m][l][f][1]->Draw("same");
        hQ[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hQ[m][l][f][0]->GetBinContent(hQ[m][l][f][0]->GetMaximumBin()), hQ[m][l][f][1]->GetBinContent(hQ[m][l][f][1]->GetMaximumBin())) + 50);
        format_h_1D(hQ[m][l][f][0]);
        format_h_1D(hQ[m][l][f][1]);
        }
//         canT->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hT[m][l][f][0]->SetLineColor(colL);
//         hT[m][l][f][1]->SetLineColor(colR);
//         hT[m][l][f][0]->Draw();
//         hT[m][l][f][1]->Draw("same");
//         hT[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hT[m][l][f][0]->GetBinContent(hT[m][l][f][0]->GetMaximumBin()), hT[m][l][f][1]->GetBinContent(hT[m][l][f][1]->GetMaximumBin())) + 50);
//         format_h_1D(hT[m][l][f][0]);
//         format_h_1D(hT[m][l][f][1]); 
//         
//         canTDif->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hTDif[m][l][f]->SetLineColor(colL);
//         hTDif[m][l][f]->Draw();
//         hTDif[m][l][f]->GetYaxis()->SetRangeUser(0,hTDif[m][l][f]->GetBinContent(hTDif[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hTDif[m][l][f]);
//         
//         canQSum->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hQSum[m][l][f]->SetLineColor(colL);
//         hQSum[m][l][f]->Draw();
//         hQSum[m][l][f]->GetYaxis()->SetRangeUser(0,hQSum[m][l][f]->GetBinContent(hQSum[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hQSum[m][l][f]);
//         
//         canQAve->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hQAve[m][l][f]->SetLineColor(colL);
//         hQAve[m][l][f]->Draw();
//         hQAve[m][l][f]->GetYaxis()->SetRangeUser(0,hQAve[m][l][f]->GetBinContent(hQAve[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hQAve[m][l][f]);
//         
//         canMLR->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hMLR[m][l][f]->SetLineColor(colL);
//         hMLR[m][l][f]->Draw();
//         hMLR[m][l][f]->GetYaxis()->SetRangeUser(0,hMLR[m][l][f]->GetBinContent(hMLR[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hMLR[m][l][f]);
//         
//         canQLvsQR->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hQLvsQR[m][l][f]->Draw("colz");
//         
//         canTLvsTR->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hTLvsTR[m][l][f]->Draw("colz");
               
//     } 

	canMult->cd();  
    gPad->SetGrid(1, 1);
	hFiberMult->Draw();
 	format_h_1D(hFiberMult);
    hFiberMult->GetXaxis()->SetRange(1, activeAddresses.size()+1);
	hFiberMult->GetXaxis()->SetTitle("# SiPMs");
	hFiberMult->GetYaxis()->SetTitle("Counts");
    
	output->mkdir("Canvases");
	output->cd("Canvases");
    
	canQDC_layer0->Write();
    canQDC_layer1->Write();
    canQDC_layer2->Write();
    canQDC_layer3->Write();
//     canT->Write();
// 	canTDif->Write();
//     canQAve->Write();
// 	canQSum->Write();
//     canMLR->Write();
// 	canQLvsQR->Write(); 
// 	canTLvsTR->Write();    
    canMult->Write();
    
//     output->Close();

// //----- writing canvas to file    
    TString pdfname = ConstructPdfName(path);
    canQDC_layer0->Print(Form("%s(",pdfname.Data()), "pdf");
    canQDC_layer1->Print(Form("%s(",pdfname.Data()), "pdf");
    canQDC_layer2->Print(Form("%s(",pdfname.Data()), "pdf");
    canQDC_layer3->Print(Form("%s(",pdfname.Data()), "pdf");
//     canT->Print(pdfname,"pdf");
//     canTDif->Print(pdfname,"pdf");
//     canQAve->Print(pdfname,"pdf");
//     canQSum->Print(pdfname,"pdf");
//     canMLR->Print(pdfname,"pdf");
//     canQLvsQR->Print(pdfname,"pdf");
//     canTLvsTR->Print(pdfname,"pdf");
    canMult->Print(Form("%s)",pdfname.Data()),"pdf");

    auto end = std::chrono::system_clock::now();   
    std::chrono::duration<double> time = end-start;
    std::cout << "Time: " << time.count() << " s\n";
    
	return 0;
}
