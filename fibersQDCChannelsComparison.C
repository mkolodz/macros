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
// #define XUP_QDC 800
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



int fibersQDCChannelsComparison(std::vector<TString> path)
{
//     std::cout << "dupksko" << std::endl;
    auto start = std::chrono::system_clock::now();  
    int ii=0;
    int fiberA=46;
    int layerA=0;
    int fiberB=54;
    int layerB=0;
    int entriesA=0;
    int entriesB=0;
    int entriesAll[10];
    for(int i=0;i<10;i++){
        entriesAll[i]=0;
    }
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

//     TH2D * hFiberHeatMap = new TH2D("hFiberHeatMap", "hFiberHeatMap", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
//     TH1D * hQDCLChanA = new TH1D("hQDCLChanA", Form("QDCL_M0L%dF%d;QDC (a.u.);Counts", layerA, fiberA), NBINS_QDC, XLOW_QDC, XUP_QDC);
//     TH1D * hQDCLChanB = new TH1D("hQDCLChanB", Form("QDCL_M0L%dF%d;QDC (a.u.);Counts", layerB, fiberB), NBINS_QDC, XLOW_QDC, XUP_QDC);
//     TH1D * hQDCRChanA = new TH1D("hQDCRChanA", Form("QDCR_M0L%dF%d;QDC (a.u.);Counts", layerA, fiberA), NBINS_QDC, XLOW_QDC, XUP_QDC);
//     TH1D * hQDCRChanB = new TH1D("hQDCRChanB", Form("QDCR_M0L%dF%d;QDC (a.u.);Counts", layerB, fiberB), NBINS_QDC, XLOW_QDC, XUP_QDC);
    TH1D * hQav[12];
    for (int k=0;k<10;k++){
       hQav[k] = new TH1D(Form("Qav_M0L%dF%d", 0, k+10), Form("Qav_M0L%dF%d;QDC (a.u.);Counts", 0, k+10), NBINS_QDC, XLOW_QDC, XUP_QDC); 
    }
    hQav[10] = new TH1D(Form("Qav_M0L%dF%d", 0, 46), Form("Qav_M0L%dF%d;QDC (a.u.);Counts", 0, 46), NBINS_QDC, XLOW_QDC, XUP_QDC); 
    hQav[11] = new TH1D(Form("Qav_M0L%dF%d", 0, 54), Form("Qav_M0L%dF%d;QDC (a.u.);Counts", 0, 54), NBINS_QDC, XLOW_QDC, XUP_QDC); 

    
//     TH1D * hQavChanB = new TH1D("hQavChanB", Form("Qav_M0L%dF%d;QDC (a.u.);Counts", layerB, fiberB), NBINS_QDC, XLOW_QDC, XUP_QDC);
	Int_t mod, lay, fi;
	SFibersRaw * pRaw;

 	int nLoop = loop->getEntries();

//----- histogram definition and filling    
	for (int i = 0; i < nLoop; ++i)
	{
		size_t nCat = pCatRaw->getEntries();
        
		for (uint j = 0; j < nCat; ++j)
		{
			pRaw = (SFibersRaw *)pCatRaw->getObject(j);
			pRaw->getAddress(mod, lay, fi);
            
            if(pRaw->getQDCL() > QDC_MIN &&
            pRaw->getQDCR() > QDC_MIN &&
            pRaw->getTimeL() > T_MIN &&
            pRaw->getTimeR() > T_MIN)
            {
//                 hFiberHeatMap->Fill(fi, lay);
                for(int kk=0;kk<10;kk++){
                    if(mod==0 && lay==0 && fi==kk+10){
                    hQav[kk]->Fill(sqrt(pRaw->getQDCR()*pRaw->getQDCL()));
                    entriesAll[kk]++;
                }
                }
                if(mod==0 && lay==layerA && fi==fiberA){
//                     hQDCLChanA->Fill(pRaw->getQDCL());
//                     hQDCRChanA->Fill(pRaw->getQDCR());
                    hQav[10]->Fill(sqrt(pRaw->getQDCR()*pRaw->getQDCL()));
                    entriesA++;
                }
                if(mod==0 && lay==layerB && fi==fiberB){
//                     hQDCLChanB->Fill(pRaw->getQDCL());
//                     hQDCRChanB->Fill(pRaw->getQDCR());
                    hQav[11]->Fill(sqrt(pRaw->getQDCR()*pRaw->getQDCL()));
                    entriesB++;
                }
            }	

		}
		
		loop->nextEvent();
	}
	
	std::cout << "\n\nLoop entries: " << nLoop << std::endl;
    TCanvas * c = new TCanvas("c", "c", 1800, 1800);
    c->Divide(3,4);
    for(int i=0;i<10;i++){
        c->cd(i+1); 
        gPad->SetGrid(1,1);
        gPad->SetLogy();
        hQav[i]->Scale(1./entriesAll[i]);
        hQav[i]->Draw("HIST");
        gPad->BuildLegend();
    }
    
    
    c->cd(11);
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    hQav[10]->Scale(1./entriesA);
    hQav[10]->Draw("HIST");
    gPad->BuildLegend();
    
    c->cd(12);
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    hQav[11]->Scale(1./entriesB);
    hQav[11]->Draw("HIST");
    gPad->BuildLegend();
    
// 	path[ii].ReplaceAll(".root", "_HEATMAPS_ONLY_TEST.root");
// 	TFile *output = new TFile(path[ii],"RECREATE");
// 	output->cd();
//     hFiberHeatMap->Write();


    auto end = std::chrono::system_clock::now();   
    std::chrono::duration<double> time = end-start;
    std::cout << "Time: " << time.count() << " s\n";
    
	return 0;
}
