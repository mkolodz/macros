// R__ADD_INCLUDE_PATH(/home/lab/sifi-framework/install/include)
R__ADD_INCLUDE_PATH(/scratch/gccb/software/framework/20230109-4to1-ubuntu20/include)
#include "SCategoryManager.h"
#include "SLoop.h"
#include "SDDSamples.h"
#include "SFibersRaw.h"
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
#define XUP_QDC 1000

#define NBINS_T 20
#define XLOW_T 0
#define XUP_T 20

#define NBINS_TDIFF 20
#define XLOW_TDIFF 0
#define XUP_TDIFF 20
// #define XUP_TDIFF_ZOOM 100

#define NBINS_MLR 100
#define XLOW_MLR -15
#define XUP_MLR 15

//----- cuts
#define TDIFF_MAX 10
#define QDC_MIN 10E-6
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

//   TString out = pathname(0,28); // length of this path: /scratch1/gccb/data/Citiroc/ - 28 characters
//   TString out = pathname(0,28); // length of this path: /scratch1/gccb/data/TOFPET2/ - 28 characters
  TString out = pathname(0,30); // 
//   TString out = pathname(0,27); // length of this path: /scratch1/gccb/data/202***/ - 27 characters
  //TString out = pathname(0,26); // length of this path: /scratch/gccb/data/202***/ - 26 characters
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




int rawDisplayFibers(std::vector<TString> path)
{
    auto start = std::chrono::system_clock::now();  
    int ii=0;
    if(!path[ii].Contains("/") || !path[ii].BeginsWith("/"))
    {
        std::cout << "##### Error! The functions needs the filename including the full absolute path..." << std::endl;
        std::abort();
    }
    
    gStyle->SetPalette(kBird);
    
    std::vector<std::string> str_path;
    for(int i=0; i < path.size(); i++){
        str_path.push_back(std::string(path[i]));
        printf("HERE: %s", str_path[i].c_str());
    }
    
	SLoop * loop = new SLoop();
	loop->addFiles(str_path);
	loop->setInput({});

	SCategory * pCatRaw = SCategoryManager::getCategory(SCategory::CatFibersRaw);
    vector<Address> activeAddresses;
    
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
    Int_t first_time_l = 0;
    Int_t first_time_r = 0;
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

            if(pRaw->getQDCL() > QDC_MIN &&
            pRaw->getQDCR() > QDC_MIN &&
            pRaw->getTimeL() > T_MIN &&
            pRaw->getTimeR() > T_MIN /*&&
            fabs(pRaw->getTimeL() - pRaw->getTimeR()) < TDIFF_MAX*/)
            {
                if(mult==0){
//                     std::cout << mult << std::endl;
                    first_time_l = pRaw->getTimeL();
                    first_time_r = pRaw->getTimeR();
                }
//                 std::cout << "tL" <<pRaw->getTimeL()-first_time_l << std::endl;
//                 std::cout  << "tR"<< pRaw->getTimeR()- first_time_r<< std::endl;
//                 std::cout  << "tDiff"<< fabs(pRaw->getTimeL() - pRaw->getTimeR() ) << std::endl;
                hQ[mod][lay][fi][0]->Fill(pRaw->getQDCL());
                hQ[mod][lay][fi][1]->Fill(pRaw->getQDCR());
                hT[mod][lay][fi][0]->Fill(pRaw->getTimeL()-first_time_l);
                hT[mod][lay][fi][1]->Fill(pRaw->getTimeR()-first_time_r);
                
                hQAve[mod][lay][fi]->Fill(sqrt(pRaw->getQDCL()*pRaw->getQDCR()));
                hQSum[mod][lay][fi]->Fill(pRaw->getQDCL() + pRaw->getQDCR());
                hMLR[mod][lay][fi]->Fill(log(sqrt(pRaw->getQDCR()/pRaw->getQDCL())));
                
                hQLvsQR[mod][lay][fi]->Fill(pRaw->getQDCL(), pRaw->getQDCR());
                hTLvsTR[mod][lay][fi]->Fill(pRaw->getTimeL()-first_time_l, pRaw->getTimeR()-first_time_r);
                hTDif[mod][lay][fi]->Fill(pRaw->getTimeL() - pRaw->getTimeR());
//                 if(fabs(pRaw->getTimeL() - pRaw->getTimeR() ) < TDIFF_MAX)
//                 {
                    
//                 }
                
            }	
++mult;
//             if(pRaw->getQDCL() > QDC_MIN &&
//             pRaw->getQDCR() > QDC_MIN &&
//             pRaw->getTimeL() > T_MIN &&
//             pRaw->getTimeR() > T_MIN)
//             {
//                 hTDif[mod][lay][fi]->Fill(pRaw->getTimeL() - pRaw->getTimeR() );
//             }
			
		}
		
		hFiberMult->Fill(mult);
		mult = 0;
		loop->nextEvent();
	}
	
	std::cout << "\n\nLoop entries: " << nLoop << std::endl;

	path[ii].ReplaceAll(".root", "_HISTOS.root");
	TFile *output = new TFile(path[ii],"RECREATE");
	output->cd();

//----- canvas definition    
//     TCanvas * canQDC_layer0 = new TCanvas("canQDC_layer0", "canQDC_layer0", XCANVAS, YCANVAS); 
//     canQDC_layer0->SetLeftMargin(0.2);
// 	canQDC_layer0->DivideSquare(activeAddresses.size()/4);
// 
//     TCanvas * canQDC_layer1 = new TCanvas("canQDC_layer1", "canQDC_layer1", XCANVAS, YCANVAS); 
//     canQDC_layer1->SetLeftMargin(0.2);
// 	canQDC_layer1->DivideSquare(activeAddresses.size()/4);
//     
//     TCanvas * canQDC_layer2 = new TCanvas("canQDC_layer2", "canQDC_layer2", XCANVAS, YCANVAS); 
//     canQDC_layer2->SetLeftMargin(0.2);
// 	canQDC_layer2->DivideSquare(activeAddresses.size()/4);
//     
//     TCanvas * canQDC_layer3 = new TCanvas("canQDC_layer3", "canQDC_layer3", XCANVAS, YCANVAS); 
//     canQDC_layer3->SetLeftMargin(0.2);
// 	canQDC_layer3->DivideSquare(activeAddresses.size()/4);
//     
//     TCanvas * canT_layer0 = new TCanvas("canT_layer0", "canT_layer0", XCANVAS, YCANVAS); 
//     canT_layer0->SetLeftMargin(0.2); 
//     canT_layer0->DivideSquare(activeAddresses.size()/4);
//     
//     TCanvas * canT_layer1 = new TCanvas("canT_layer1", "canT_layer1", XCANVAS, YCANVAS); 
//     canT_layer1->SetLeftMargin(0.2); 
//     canT_layer1->DivideSquare(activeAddresses.size()/4);
//     
//     TCanvas * canT_layer2 = new TCanvas("canT_layer2", "canT_layer2", XCANVAS, YCANVAS); 
//     canT_layer2->SetLeftMargin(0.2); 
//     canT_layer2->DivideSquare(activeAddresses.size()/4);
//     
//     TCanvas * canT_layer3 = new TCanvas("canT_layer3", "canT_layer3", XCANVAS, YCANVAS); 
//     canT_layer3->SetLeftMargin(0.2); 
//     canT_layer3->DivideSquare(activeAddresses.size()/4);
//   
// 	TCanvas * canTDif_layer0 = new TCanvas("canTDif_layer0", "canTDif_layer0", XCANVAS, YCANVAS);
//     canTDif_layer0->SetLeftMargin(0.2);
// 	canTDif_layer0->DivideSquare(activeAddresses.size()/4);
//     
// 	TCanvas * canTDif_layer1 = new TCanvas("canTDif_layer1", "canTDif_layer1", XCANVAS, YCANVAS);
//     canTDif_layer1->SetLeftMargin(0.2);
// 	canTDif_layer1->DivideSquare(activeAddresses.size()/4);
//     
// 	TCanvas * canTDif_layer2 = new TCanvas("canTDif_layer2", "canTDif_layer2", XCANVAS, YCANVAS);
//     canTDif_layer2->SetLeftMargin(0.2);
// 	canTDif_layer2->DivideSquare(activeAddresses.size()/4);
//     
// 	TCanvas * canTDif_layer3 = new TCanvas("canTDif_layer3", "canTDif_layer3", XCANVAS, YCANVAS);
//     canTDif_layer3->SetLeftMargin(0.2);
// 	canTDif_layer3->DivideSquare(activeAddresses.size()/4);
//     
// 
//     TCanvas * canQAve_layer0 = new TCanvas("canQAve_layer0", "canQAve_layer0", XCANVAS, YCANVAS); 
//     canQAve_layer0->SetLeftMargin(0.2); 
//     canQAve_layer0->DivideSquare(activeAddresses.size()/4);
//     
//     TCanvas * canQAve_layer1 = new TCanvas("canQAve_layer1", "canQAve_layer1", XCANVAS, YCANVAS); 
//     canQAve_layer1->SetLeftMargin(0.2); 
//     canQAve_layer1->DivideSquare(activeAddresses.size()/4);
//     
//     TCanvas * canQAve_layer2 = new TCanvas("canQAve_layer2", "canQAve_layer2", XCANVAS, YCANVAS); 
//     canQAve_layer2->SetLeftMargin(0.2); 
//     canQAve_layer2->DivideSquare(activeAddresses.size()/4);
//     
//     TCanvas * canQAve_layer3 = new TCanvas("canQAve_layer3", "canQAve_layer3", XCANVAS, YCANVAS); 
//     canQAve_layer3->SetLeftMargin(0.2); 
//     canQAve_layer3->DivideSquare(activeAddresses.size()/4);
//     
//     
// 	TCanvas * canQSum_layer0 = new TCanvas("canQSum_layer0", "canQSum_layer0", XCANVAS, YCANVAS); 
//     canQSum_layer0->SetLeftMargin(0.2); 
// 	canQSum_layer0->DivideSquare(activeAddresses.size()/4);
//     
//     	TCanvas * canQSum_layer1 = new TCanvas("canQSum_layer1", "canQSum_layer1", XCANVAS, YCANVAS); 
//     canQSum_layer1->SetLeftMargin(0.2); 
// 	canQSum_layer1->DivideSquare(activeAddresses.size()/4);
//     
//     	TCanvas * canQSum_layer2 = new TCanvas("canQSum_layer2", "canQSum_layer2", XCANVAS, YCANVAS); 
//     canQSum_layer2->SetLeftMargin(0.2); 
// 	canQSum_layer2->DivideSquare(activeAddresses.size()/4);
//     
//     	TCanvas * canQSum_layer3 = new TCanvas("canQSum_layer3", "canQSum_layer3", XCANVAS, YCANVAS); 
//     canQSum_layer3->SetLeftMargin(0.2); 
// 	canQSum_layer3->DivideSquare(activeAddresses.size()/4);
//     
// 	TCanvas * canMLR_layer0 = new TCanvas("canMLR_layer0", "canMLR_layer0", XCANVAS, YCANVAS); 
//     canMLR_layer0->SetLeftMargin(0.2); 
// 	canMLR_layer0->DivideSquare(activeAddresses.size()/4);
//         
// 	TCanvas * canMLR_layer1 = new TCanvas("canMLR_layer1", "canMLR_layer1", XCANVAS, YCANVAS); 
//     canMLR_layer1->SetLeftMargin(0.2); 
// 	canMLR_layer1->DivideSquare(activeAddresses.size()/4);
//         
// 	TCanvas * canMLR_layer2 = new TCanvas("canMLR_layer2", "canMLR_layer2", XCANVAS, YCANVAS); 
//     canMLR_layer2->SetLeftMargin(0.2); 
// 	canMLR_layer2->DivideSquare(activeAddresses.size()/4);
//         
// 	TCanvas * canMLR_layer3 = new TCanvas("canMLR_layer3", "canMLR_layer3", XCANVAS, YCANVAS); 
//     canMLR_layer3->SetLeftMargin(0.2); 
// 	canMLR_layer3->DivideSquare(activeAddresses.size()/4);
//     
//     
// 	TCanvas * canQLvsQR_layer0 = new TCanvas("canQLvsQR_layer0", "canQLvsQR_layer0", XCANVAS, YCANVAS); 
//     canQLvsQR_layer0->SetLeftMargin(0.2);
// 	canQLvsQR_layer0->DivideSquare(activeAddresses.size()/4);
//         
// 	TCanvas * canQLvsQR_layer1 = new TCanvas("canQLvsQR_layer1", "canQLvsQR_layer1", XCANVAS, YCANVAS); 
//     canQLvsQR_layer1->SetLeftMargin(0.2);
// 	canQLvsQR_layer1->DivideSquare(activeAddresses.size()/4);
//     
//     
// 	TCanvas * canQLvsQR_layer2 = new TCanvas("canQLvsQR_layer2", "canQLvsQR_layer2", XCANVAS, YCANVAS); 
//     canQLvsQR_layer2->SetLeftMargin(0.2);
// 	canQLvsQR_layer2->DivideSquare(activeAddresses.size()/4);
//         
// 	TCanvas * canQLvsQR_layer3 = new TCanvas("canQLvsQR_layer3", "canQLvsQR_layer3", XCANVAS, YCANVAS); 
//     canQLvsQR_layer3->SetLeftMargin(0.2);
// 	canQLvsQR_layer3->DivideSquare(activeAddresses.size()/4);
//     
//     
//     
//     TCanvas * canTLvsTR_layer0 = new TCanvas("canTLvsTR_layer0", "canTLvsTR_layer0", XCANVAS, YCANVAS); 
//     canTLvsTR_layer0->SetLeftMargin(0.2);
// 	canTLvsTR_layer0->DivideSquare(activeAddresses.size()/4);
//     
//     
//     TCanvas * canTLvsTR_layer1 = new TCanvas("canTLvsTR_layer1", "canTLvsTR_layer1", XCANVAS, YCANVAS); 
//     canTLvsTR_layer1->SetLeftMargin(0.2);
// 	canTLvsTR_layer1->DivideSquare(activeAddresses.size()/4);
//         
//     
//     TCanvas * canTLvsTR_layer2 = new TCanvas("canTLvsTR_layer2", "canTLvsTR_layer2", XCANVAS, YCANVAS); 
//     canTLvsTR_layer2->SetLeftMargin(0.2);
// 	canTLvsTR_layer2->DivideSquare(activeAddresses.size()/4);
//         
//     
//     TCanvas * canTLvsTR_layer3 = new TCanvas("canTLvsTR_layer3", "canTLvsTR_layer3", XCANVAS, YCANVAS); 
//     canTLvsTR_layer3->SetLeftMargin(0.2);
// 	canTLvsTR_layer3->DivideSquare(activeAddresses.size()/4);
//     
// 	TCanvas * canMult = new TCanvas("Mult", "Mult", 600, 400); 
//     canMult->SetLeftMargin(0.2);
    
    
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

        hTDif[m][l][f]->SetDirectory(gDirectory);
        hTDif[m][l][f]->Write();
        hQSum[m][l][f]->SetDirectory(gDirectory);
        hQSum[m][l][f]->Write();
        hQAve[m][l][f]->SetDirectory(gDirectory);
        hQAve[m][l][f]->Write();
        hMLR[m][l][f]->SetDirectory(gDirectory);
        hMLR[m][l][f]->Write();
        
        hQLvsQR[m][l][f]->SetDirectory(gDirectory);
        hQLvsQR[m][l][f]->Write();
        hTLvsTR[m][l][f]->SetDirectory(gDirectory);
        hTLvsTR[m][l][f]->Write();
        

//         canQDC_layer0->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hQ[m][l][f][0]->SetLineColor(colL);
//         hQ[m][l][f][1]->SetLineColor(colR);
//         hQ[m][l][f][0]->Draw();
//         hQ[m][l][f][1]->Draw("same");
//         hQ[m][l][f][0]->GetYaxis()->SetRangeUser(0.1, TMath::Max(hQ[m][l][f][0]->GetBinContent(hQ[m][l][f][0]->GetMaximumBin()), hQ[m][l][f][1]->GetBinContent(hQ[m][l][f][1]->GetMaximumBin())) + 50);
//         gPad->SetLogy();
//         format_h_1D(hQ[m][l][f][0]);
//         format_h_1D(hQ[m][l][f][1]);
//         
//         canT_layer0->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hT[m][l][f][0]->SetLineColor(colL);
//         hT[m][l][f][1]->SetLineColor(colR);
//         hT[m][l][f][0]->Draw();
//         hT[m][l][f][1]->Draw("same");
//         hT[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hT[m][l][f][0]->GetBinContent(hT[m][l][f][0]->GetMaximumBin()), hT[m][l][f][1]->GetBinContent(hT[m][l][f][1]->GetMaximumBin())) + 50);
//         format_h_1D(hT[m][l][f][0]);
//         format_h_1D(hT[m][l][f][1]); 
// 
//         canTDif_layer0->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hTDif[m][l][f]->SetLineColor(colL);
//         hTDif[m][l][f]->Draw();
//         hTDif[m][l][f]->GetYaxis()->SetRangeUser(0,hTDif[m][l][f]->GetBinContent(hTDif[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hTDif[m][l][f]);
//         
//         canQSum_layer0->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hQSum[m][l][f]->SetLineColor(colL);
//         hQSum[m][l][f]->Draw();
//         hQSum[m][l][f]->GetYaxis()->SetRangeUser(0,hQSum[m][l][f]->GetBinContent(hQSum[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hQSum[m][l][f]);
//         
//         canQAve_layer0->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hQAve[m][l][f]->SetLineColor(colL);
//         hQAve[m][l][f]->Draw();
//         hQAve[m][l][f]->GetYaxis()->SetRangeUser(0,hQAve[m][l][f]->GetBinContent(hQAve[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hQAve[m][l][f]);
//         
//         canMLR_layer0->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hMLR[m][l][f]->SetLineColor(colL);
//         hMLR[m][l][f]->Draw();
//         hMLR[m][l][f]->GetYaxis()->SetRangeUser(0,hMLR[m][l][f]->GetBinContent(hMLR[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hMLR[m][l][f]);
//         
//         canQLvsQR_layer0->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hQLvsQR[m][l][f]->Draw("colz");
//         
//         canTLvsTR_layer0->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hTLvsTR[m][l][f]->Draw("colz");
        
        
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

        hTDif[m][l][f]->SetDirectory(gDirectory);
        hTDif[m][l][f]->Write();
        hQSum[m][l][f]->SetDirectory(gDirectory);
        hQSum[m][l][f]->Write();
        hQAve[m][l][f]->SetDirectory(gDirectory);
        hQAve[m][l][f]->Write();
        hMLR[m][l][f]->SetDirectory(gDirectory);
        hMLR[m][l][f]->Write();
        
        hQLvsQR[m][l][f]->SetDirectory(gDirectory);
        hQLvsQR[m][l][f]->Write();
        hTLvsTR[m][l][f]->SetDirectory(gDirectory);
        hTLvsTR[m][l][f]->Write();
        
//         canQDC_layer1->cd(i-activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hQ[m][l][f][0]->SetLineColor(colL);
//         hQ[m][l][f][1]->SetLineColor(colR);
//         hQ[m][l][f][0]->Draw();
//         hQ[m][l][f][1]->Draw("same");
//         hQ[m][l][f][0]->GetYaxis()->SetRangeUser(0.1, TMath::Max(hQ[m][l][f][0]->GetBinContent(hQ[m][l][f][0]->GetMaximumBin()), hQ[m][l][f][1]->GetBinContent(hQ[m][l][f][1]->GetMaximumBin())) + 50);
//         gPad->SetLogy();
//         format_h_1D(hQ[m][l][f][0]);
//         format_h_1D(hQ[m][l][f][1]);
//         
//         canT_layer1->cd(i-activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hT[m][l][f][0]->SetLineColor(colL);
//         hT[m][l][f][1]->SetLineColor(colR);
//         hT[m][l][f][0]->Draw();
//         hT[m][l][f][1]->Draw("same");
//         hT[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hT[m][l][f][0]->GetBinContent(hT[m][l][f][0]->GetMaximumBin()), hT[m][l][f][1]->GetBinContent(hT[m][l][f][1]->GetMaximumBin())) + 50);
//         format_h_1D(hT[m][l][f][0]);
//         format_h_1D(hT[m][l][f][1]); 
// 
//         canTDif_layer1->cd(i-activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hTDif[m][l][f]->SetLineColor(colL);
//         hTDif[m][l][f]->Draw();
//         hTDif[m][l][f]->GetYaxis()->SetRangeUser(0,hTDif[m][l][f]->GetBinContent(hTDif[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hTDif[m][l][f]);
//         
//         canQSum_layer1->cd(i-activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hQSum[m][l][f]->SetLineColor(colL);
//         hQSum[m][l][f]->Draw();
//         hQSum[m][l][f]->GetYaxis()->SetRangeUser(0,hQSum[m][l][f]->GetBinContent(hQSum[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hQSum[m][l][f]);
//         
//         canQAve_layer1->cd(i-activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hQAve[m][l][f]->SetLineColor(colL);
//         hQAve[m][l][f]->Draw();
//         hQAve[m][l][f]->GetYaxis()->SetRangeUser(0,hQAve[m][l][f]->GetBinContent(hQAve[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hQAve[m][l][f]);
//         
//         canMLR_layer1->cd(i-activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hMLR[m][l][f]->SetLineColor(colL);
//         hMLR[m][l][f]->Draw();
//         hMLR[m][l][f]->GetYaxis()->SetRangeUser(0,hMLR[m][l][f]->GetBinContent(hMLR[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hMLR[m][l][f]);
//         
//         canQLvsQR_layer1->cd(i-activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hQLvsQR[m][l][f]->Draw("colz");
//         
//         canTLvsTR_layer1->cd(i-activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hTLvsTR[m][l][f]->Draw("colz");
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

        hTDif[m][l][f]->SetDirectory(gDirectory);
        hTDif[m][l][f]->Write();
        hQSum[m][l][f]->SetDirectory(gDirectory);
        hQSum[m][l][f]->Write();
        hQAve[m][l][f]->SetDirectory(gDirectory);
        hQAve[m][l][f]->Write();
        hMLR[m][l][f]->SetDirectory(gDirectory);
        hMLR[m][l][f]->Write();
        
        hQLvsQR[m][l][f]->SetDirectory(gDirectory);
        hQLvsQR[m][l][f]->Write();
        hTLvsTR[m][l][f]->SetDirectory(gDirectory);
        hTLvsTR[m][l][f]->Write();
//         canQDC_layer2->cd(i-activeAddresses.size()/2+1);
//         gPad->SetGrid(1, 1);
//         hQ[m][l][f][0]->SetLineColor(colL);
//         hQ[m][l][f][1]->SetLineColor(colR);
//         hQ[m][l][f][0]->Draw();
//         hQ[m][l][f][1]->Draw("same");
//         hQ[m][l][f][0]->GetYaxis()->SetRangeUser(0.1, TMath::Max(hQ[m][l][f][0]->GetBinContent(hQ[m][l][f][0]->GetMaximumBin()), hQ[m][l][f][1]->GetBinContent(hQ[m][l][f][1]->GetMaximumBin())) + 50);        
//         gPad->SetLogy();
//         format_h_1D(hQ[m][l][f][0]);
//         format_h_1D(hQ[m][l][f][1]);
//         
//         canT_layer2->cd(i-activeAddresses.size()/2+1);
//         gPad->SetGrid(1, 1);
//         hT[m][l][f][0]->SetLineColor(colL);
//         hT[m][l][f][1]->SetLineColor(colR);
//         hT[m][l][f][0]->Draw();
//         hT[m][l][f][1]->Draw("same");
//         hT[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hT[m][l][f][0]->GetBinContent(hT[m][l][f][0]->GetMaximumBin()), hT[m][l][f][1]->GetBinContent(hT[m][l][f][1]->GetMaximumBin())) + 50);
//         format_h_1D(hT[m][l][f][0]);
//         format_h_1D(hT[m][l][f][1]); 
// 
//         canTDif_layer2->cd(i-activeAddresses.size()/2+1);
//         gPad->SetGrid(1, 1);
//         hTDif[m][l][f]->SetLineColor(colL);
//         hTDif[m][l][f]->Draw();
//         hTDif[m][l][f]->GetYaxis()->SetRangeUser(0,hTDif[m][l][f]->GetBinContent(hTDif[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hTDif[m][l][f]);
//         
//         canQSum_layer2->cd(i-activeAddresses.size()/2+1);
//         gPad->SetGrid(1, 1);
//         hQSum[m][l][f]->SetLineColor(colL);
//         hQSum[m][l][f]->Draw();
//         hQSum[m][l][f]->GetYaxis()->SetRangeUser(0,hQSum[m][l][f]->GetBinContent(hQSum[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hQSum[m][l][f]);
//         
//         canQAve_layer2->cd(i-activeAddresses.size()/2+1);
//         gPad->SetGrid(1, 1);
//         hQAve[m][l][f]->SetLineColor(colL);
//         hQAve[m][l][f]->Draw();
//         hQAve[m][l][f]->GetYaxis()->SetRangeUser(0,hQAve[m][l][f]->GetBinContent(hQAve[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hQAve[m][l][f]);
//         
//         canMLR_layer2->cd(i-activeAddresses.size()/2+1);
//         gPad->SetGrid(1, 1);
//         hMLR[m][l][f]->SetLineColor(colL);
//         hMLR[m][l][f]->Draw();
//         hMLR[m][l][f]->GetYaxis()->SetRangeUser(0,hMLR[m][l][f]->GetBinContent(hMLR[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hMLR[m][l][f]);
//         
//         canQLvsQR_layer2->cd(i-activeAddresses.size()/2+1);
//         gPad->SetGrid(1, 1);
//         hQLvsQR[m][l][f]->Draw("colz");
//         
//         canTLvsTR_layer2->cd(i-activeAddresses.size()/2+1);
//         gPad->SetGrid(1, 1);
//         hTLvsTR[m][l][f]->Draw("colz");
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

        hTDif[m][l][f]->SetDirectory(gDirectory);
        hTDif[m][l][f]->Write();
        hQSum[m][l][f]->SetDirectory(gDirectory);
        hQSum[m][l][f]->Write();
        hQAve[m][l][f]->SetDirectory(gDirectory);
        hQAve[m][l][f]->Write();
        hMLR[m][l][f]->SetDirectory(gDirectory);
        hMLR[m][l][f]->Write();
        
        hQLvsQR[m][l][f]->SetDirectory(gDirectory);
        hQLvsQR[m][l][f]->Write();
        hTLvsTR[m][l][f]->SetDirectory(gDirectory);
        hTLvsTR[m][l][f]->Write();
        
//         canQDC_layer3->cd(i-3*activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hQ[m][l][f][0]->SetLineColor(colL);
//         hQ[m][l][f][1]->SetLineColor(colR);
//         hQ[m][l][f][0]->Draw();
//         hQ[m][l][f][1]->Draw("same");
//         hQ[m][l][f][0]->GetYaxis()->SetRangeUser(0.1, TMath::Max(hQ[m][l][f][0]->GetBinContent(hQ[m][l][f][0]->GetMaximumBin()), hQ[m][l][f][1]->GetBinContent(hQ[m][l][f][1]->GetMaximumBin())) + 50);
//         gPad->SetLogy();
//         format_h_1D(hQ[m][l][f][0]);
//         format_h_1D(hQ[m][l][f][1]);
//         
//         canT_layer3->cd(i-3*activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hT[m][l][f][0]->SetLineColor(colL);
//         hT[m][l][f][1]->SetLineColor(colR);
//         hT[m][l][f][0]->Draw();
//         hT[m][l][f][1]->Draw("same");
//         hT[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hT[m][l][f][0]->GetBinContent(hT[m][l][f][0]->GetMaximumBin()), hT[m][l][f][1]->GetBinContent(hT[m][l][f][1]->GetMaximumBin())) + 50);
//         format_h_1D(hT[m][l][f][0]);
//         format_h_1D(hT[m][l][f][1]); 
// 
//         canTDif_layer3->cd(i-3*activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hTDif[m][l][f]->SetLineColor(colL);
//         hTDif[m][l][f]->Draw();
//         hTDif[m][l][f]->GetYaxis()->SetRangeUser(0,hTDif[m][l][f]->GetBinContent(hTDif[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hTDif[m][l][f]);
//         
//         canQSum_layer3->cd(i-3*activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hQSum[m][l][f]->SetLineColor(colL);
//         hQSum[m][l][f]->Draw();
//         hQSum[m][l][f]->GetYaxis()->SetRangeUser(0,hQSum[m][l][f]->GetBinContent(hQSum[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hQSum[m][l][f]);
//         
//         canQAve_layer3->cd(i-3*activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hQAve[m][l][f]->SetLineColor(colL);
//         hQAve[m][l][f]->Draw();
//         hQAve[m][l][f]->GetYaxis()->SetRangeUser(0,hQAve[m][l][f]->GetBinContent(hQAve[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hQAve[m][l][f]);
//         
//         canMLR_layer3->cd(i-3*activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hMLR[m][l][f]->SetLineColor(colL);
//         hMLR[m][l][f]->Draw();
//         hMLR[m][l][f]->GetYaxis()->SetRangeUser(0,hMLR[m][l][f]->GetBinContent(hMLR[m][l][f]->GetMaximumBin()) + 50);
//         format_h_1D(hMLR[m][l][f]);
//         
//         canQLvsQR_layer3->cd(i-3*activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hQLvsQR[m][l][f]->Draw("colz");
//         
//         canTLvsTR_layer3->cd(i-3*activeAddresses.size()/4+1);
//         gPad->SetGrid(1, 1);
//         hTLvsTR[m][l][f]->Draw("colz");
        }


// 	canMult->cd();  
//     gPad->SetGrid(1, 1);
// 	hFiberMult->Draw();
//  	format_h_1D(hFiberMult);
//     hFiberMult->GetXaxis()->SetRange(1, activeAddresses.size()+1);
// 	hFiberMult->GetXaxis()->SetTitle("# SiPMs");
// 	hFiberMult->GetYaxis()->SetTitle("Counts");
//     
// 	output->mkdir("Canvases");
// 	output->cd("Canvases");
    
// 	canQDC_layer0->Write();
//     canQDC_layer1->Write();
//     canQDC_layer2->Write();
//     canQDC_layer3->Write();
//     canT_layer0->Write();
//     canT_layer1->Write();
//     canT_layer2->Write();
//     canT_layer3->Write();
// 	canTDif_layer0->Write();
// 	canTDif_layer1->Write();
// 	canTDif_layer2->Write();
// 	canTDif_layer3->Write();
//     canQAve_layer0->Write();
//     canQAve_layer1->Write();
//     canQAve_layer2->Write();
//     canQAve_layer3->Write();
// 	canQSum_layer0->Write();
// 	canQSum_layer1->Write();
// 	canQSum_layer2->Write();
// 	canQSum_layer3->Write();
//     canMLR_layer0->Write();
//     canMLR_layer1->Write();
//     canMLR_layer2->Write();
//     canMLR_layer3->Write();
// 	canQLvsQR_layer0->Write(); 
// 	canQLvsQR_layer1->Write(); 
// 	canQLvsQR_layer2->Write(); 
// 	canQLvsQR_layer3->Write(); 
// 	canTLvsTR_layer0->Write();    
// 	canTLvsTR_layer1->Write();    
// 	canTLvsTR_layer2->Write();    
// 	canTLvsTR_layer3->Write();    
//     canMult->Write();
    
//     output->Close();

// //----- writing canvas to file    
//     TString pdfname = ConstructPdfName(path[ii]);
//     canQDC_layer0->Print(Form("%s(",pdfname.Data()), "pdf");
//     canQDC_layer1->Print(Form("%s(",pdfname.Data()), "pdf");
//     canQDC_layer2->Print(Form("%s(",pdfname.Data()), "pdf");
//     canQDC_layer3->Print(Form("%s(",pdfname.Data()), "pdf");
//     canT_layer0->Print(pdfname,"pdf");
//     canT_layer1->Print(pdfname,"pdf");
//     canT_layer2->Print(pdfname,"pdf");
//     canT_layer3->Print(pdfname,"pdf");
//     canTDif_layer0->Print(pdfname,"pdf");
//     canTDif_layer1->Print(pdfname,"pdf");
//     canTDif_layer2->Print(pdfname,"pdf");
//     canTDif_layer3->Print(pdfname,"pdf");
//     canQAve_layer0->Print(pdfname,"pdf");
//     canQAve_layer1->Print(pdfname,"pdf");
//     canQAve_layer2->Print(pdfname,"pdf");
//     canQAve_layer3->Print(pdfname,"pdf");
//     canQSum_layer0->Print(pdfname,"pdf");
//     canQSum_layer1->Print(pdfname,"pdf");
//     canQSum_layer2->Print(pdfname,"pdf");
//     canQSum_layer3->Print(pdfname,"pdf");
//     canMLR_layer0->Print(pdfname,"pdf");
//     canMLR_layer1->Print(pdfname,"pdf");
//     canMLR_layer2->Print(pdfname,"pdf");
//     canMLR_layer3->Print(pdfname,"pdf");
//     canQLvsQR_layer0->Print(pdfname,"pdf");
//     canQLvsQR_layer1->Print(pdfname,"pdf");
//     canQLvsQR_layer2->Print(pdfname,"pdf");
//     canQLvsQR_layer3->Print(pdfname,"pdf");
//     canTLvsTR_layer0->Print(pdfname,"pdf");
//     canTLvsTR_layer1->Print(pdfname,"pdf");
//     canTLvsTR_layer2->Print(pdfname,"pdf");
//     canTLvsTR_layer3->Print(pdfname,"pdf");
//     canMult->Print(Form("%s)",pdfname.Data()),"pdf");

    auto end = std::chrono::system_clock::now();   
    std::chrono::duration<double> time = end-start;
    std::cout << "Time: " << time.count() << " s\n";
    
	return 0;
}
