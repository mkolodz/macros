#include "SCategoryManager.h"
#include "SLoop.h"
#include "SDDSamples.h"
#include "SFibersRaw.h"
// #include "SiFiSeries.h"
// #include "SiFiDatabase.h"
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
// #include <cstring>

//----- prototype geometry
#define N_MODULES 2
#define N_LAYERS_PER_MODULE 30
#define N_FIBERS_PER_LAYER 76
#define N_SIDES 2

//----- binning and histogram ranges
#define NBINS_ADC 200
#define XLOW_ADC 0
#define XUP_ADC 100

#define NBINS_QDC 250
#define XLOW_QDC 0
#define XUP_QDC 3000

#define NBINS_T 500
#define XLOW_T 0
#define XUP_T 500

#define NBINS_TDIFF 200
#define XLOW_TDIFF -50
#define XUP_TDIFF 50

#define NBINS_MLR 400
#define XLOW_MLR -5
#define XUP_MLR 5

//----- cuts
#define ADC_MAX 600
#define TDIFF_MAX 10
#define QDC_MIN 0
#define TOT_MIN 0
#define T_MIN 0

//----- canvas size
#define XCANVAS 1600
#define YCANVAS 1000

//----- auxiliary - to simplify the analysis when we have known # of active channels (it is in the data folder, in the params.txt/[FibersDDLookupTable]), works only if this is a continuous sequence of fibers (e.g. fibers 0-5 from the given M and L, not fibers 0-3 and then 10-15). Can be necessary to change this when changing daya series.
#define FIRST_ACTIVE_FIBER 0
#define LAST_ACTIVE_FIBER 7
#define FIXED_MODULE 0
#define FIXED_LAYER 0
#define N_ACTIVE_FIBERS 8

/*! \file rawDisplayDynamic.C
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

  TString out = pathname(0,27); // length of this path: /scratch1/gccb/data/202***/ - 27 characters
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

void PaintBin (TH1 *h, Int_t bin, Int_t color) 
{
//    printf("%d %d %d\n", bin, color, h->GetBinContent(bin));
   TBox *b = new TBox(h->GetBinLowEdge(bin),
                      h->GetMinimum(),
                      h->GetBinWidth(bin)+h->GetBinLowEdge(bin),
                      h->GetBinContent(bin));
   b->SetFillStyle(3002);
   b->SetFillColor(38);
   b->Draw();

}

int rawDisplayDynamicCrosstalk(TString path)
{
    auto start = std::chrono::system_clock::now();  
    
    if(!path.Contains("/") || !path.BeginsWith("/"))
    {
        std::cout << "##### Error! The functions needs the filename including the full absolute path..." << std::endl;
        std::abort();
    }
    
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);
	SLoop * loop = new SLoop();
	loop->addFile(std::string(path));
	loop->setInput({});

	SCategory * pCatRaw = SCategoryManager::getCategory(SCategory::CatFibersRaw);
	SCategory * pCatSample = SCategoryManager::getCategory(SCategory::CatDDSamples);
    
    struct Address {
        int iMod, iLay, iFib;
    }; 
    
    vector<Address> activeAddresses;
    int auxiliaryCounts[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES] = {}; //change to bool?
    double auxiliaryQvalues[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES]= {};
    
//     SiFiDatabase * db = new SiFiDatabase("/path/to/database.db");
//     SiFiSeries * sifi = db->GetSeries(SERIES_NO)
//     int N_MODULES = sifi->fModules;
//     for(int i = 0; i< N_MODULES; i++){
// //         sifi->fLayers;    
//     }
    
    TH1D * hA[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES];
	TH1D * hQ[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES];
    TH1D * hT[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES];
    TH1D * hCrosstalk[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES];
    
	TH1D * hQSum[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER];
    TH1D * hQAve[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER];
	TH1D * hTDif[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER];
    TH1D * hMLR[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER];
    
	TH2D * hQLvsQR[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER];
    TH2D * hTLvsTR[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER];
    
	TH1D * hFiberMult = new TH1D("hFiberMult", "hFiberMult", N_FIBERS_PER_LAYER+1, -0.5, N_FIBERS_PER_LAYER+0.5);
    TH2D * hQn1NormvsQn[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES];
    TH2D * hQn_1NormvsQn[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES];


	Int_t mod, lay, fi, side;
	SDDSamples * samples;
	SFibersRaw * pRaw;
	SDDSignal  * sigL;
	SDDSignal  * sigR;
	UShort_t mult = 0;
    int nCoincidences = 0;
    int maxDepositionFiber = 0;
    Bool_t gapBetweenNeighbours = 0;
    
 	int nLoop = loop->getEntries();

//----- histogram definition and filling    
	for (int i = 0; i < nLoop; ++i)
	{
		size_t nCat = pCatRaw->getEntries();
        
		for (uint j = 0; j < nCat; ++j)
		{
			samples = (SDDSamples*)pCatSample->getObject(j);
			pRaw = (SFibersRaw *)pCatRaw->getObject(j);
			sigL = (SDDSignal*)samples->getSignalL();
			sigR = (SDDSignal*)samples->getSignalR();
			pRaw->getAddress(mod, lay, fi);
            
            if(!hQ[mod][lay][fi][0] /*|| !hA[mod][lay][fi][0]*/){
                std::string suffix = std::string("_M") + std::to_string(mod) + std::string("L") + std::to_string(lay) + std::string("F") + std::to_string(fi);
                std::string A_L = std::string("A") + suffix + std::string("L");
                std::string A_R = std::string("A") + suffix + std::string("R");
				std::string Q_L = std::string("Q") + suffix + std::string("L");
				std::string Q_R = std::string("Q") + suffix + std::string("R");
				std::string T_L = std::string("T") + suffix + std::string("L");
				std::string T_R = std::string("T") + suffix + std::string("R");
                std::string Crosstalk_L = std::string("Crosstalk") + suffix + std::string("L");
				std::string Crosstalk_R = std::string("Crosstalk") + suffix + std::string("R");
				std::string TDif = std::string("TDif") + suffix;
                std::string QAve = std::string("QAve") + suffix;
				std::string QSum = std::string("QSum") + suffix;
				std::string MLR = std::string("MLR") + suffix;
				std::string QLvsQR = std::string("QLvsQR") + suffix;
                std::string TLvsTR = std::string("TLvsTR") + suffix;
                std::string Qn1NormvsQn_L = std::string("Qn1NormvsQnL") + suffix;
                std::string Qn1NormvsQn_R = std::string("Qn1NormvsQnR") + suffix;
                std::string Qn_1NormvsQn_L = std::string("Qn_1NormvsQnL") + suffix;
                std::string Qn_1NormvsQn_R = std::string("Qn_1NormvsQnR") + suffix;

				std::string title_A_L = A_L + std::string("; A_{L} [mV]; Counts");
				std::string title_A_R = A_R + std::string("; A_{R} [mV]; Counts");
				std::string title_Q_L = Q_L + std::string("; Q_{L} [a.u.]; Counts");
				std::string title_Q_R = Q_R + std::string("; Q_{R} [a.u.]; Counts");
				std::string title_T_L = T_L + std::string("; T_{L} [ns]; Counts");
				std::string title_T_R = T_R + std::string("; T_{R} [ns]; Counts");
                std::string title_Crosstalk_L = Crosstalk_L + std::string("; Fiber number; Counts");
				std::string title_Crosstalk_R = Crosstalk_R + std::string("; Fiber number; Counts");
				std::string title_TDif = TDif + std::string("; T_{dif} [ns]; Counts");
                std::string title_QAve = QAve + std::string("; Q_{ave} [a.u.]; Counts");
				std::string title_QSum = QSum + std::string("; Q_{L} + Q_{R} [a.u.]; Counts");
                std::string title_MLR = MLR + std::string("; M_{LR}; Counts");
				std::string title_QLvsQR = QLvsQR + std::string("; Q_{L} [a.u.]; Q_{R} [a.u.]");
				std::string title_TLvsTR = TLvsTR + std::string("; T_{L} [ns]; T_{R} [ns]");
				std::string title_Qn1NormvsQn_L = Qn1NormvsQn_L + std::string("; Q_{n+1}/Q_{n}; Q_{n} [a.u.]");
				std::string title_Qn1NormvsQn_R = Qn1NormvsQn_R + std::string("; Q_{n+1}/Q_{n}; Q_{n} [a.u.]");
				std::string title_Qn_1NormvsQn_L = Qn_1NormvsQn_L + std::string("; Q_{n-1}/Q_{n}; Q_{n} [a.u.]");
				std::string title_Qn_1NormvsQn_R = Qn_1NormvsQn_R + std::string("; Q_{n-1}/Q_{n}; Q_{n} [a.u.]");

				hA[mod][lay][fi][0] = new TH1D(A_L.c_str(), title_A_L.c_str(), NBINS_ADC, XLOW_ADC, XUP_ADC);
				hA[mod][lay][fi][1] = new TH1D(A_R.c_str(), title_A_R.c_str(), NBINS_ADC, XLOW_ADC, XUP_ADC);
				hQ[mod][lay][fi][0] = new TH1D(Q_L.c_str(), title_Q_L.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);
				hQ[mod][lay][fi][1] = new TH1D(Q_R.c_str(), title_Q_R.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);
				hT[mod][lay][fi][0] = new TH1D(T_L.c_str(), title_T_L.c_str(), NBINS_T, XLOW_T, XUP_T);
				hT[mod][lay][fi][1] = new TH1D(T_R.c_str(), title_T_R.c_str(), NBINS_T, XLOW_T, XUP_T);
                hCrosstalk[mod][lay][fi][0] = new TH1D(Crosstalk_L.c_str(), title_Crosstalk_L.c_str(), N_ACTIVE_FIBERS+1, -0.5, N_ACTIVE_FIBERS+0.5); 
                hCrosstalk[mod][lay][fi][1] = new TH1D(Crosstalk_R.c_str(), title_Crosstalk_R.c_str(), N_ACTIVE_FIBERS+1, -0.5, N_ACTIVE_FIBERS+0.5); 
				
                hTDif[mod][lay][fi] = new TH1D(TDif.c_str(), title_TDif.c_str(), NBINS_TDIFF, XLOW_TDIFF, XUP_TDIFF);
                hQAve[mod][lay][fi] = new TH1D(QAve.c_str(), title_QAve.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);
				hQSum[mod][lay][fi] = new TH1D(QSum.c_str(), title_QSum.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC*2);
                hMLR[mod][lay][fi] = new TH1D(MLR.c_str(), title_MLR.c_str(), NBINS_MLR, XLOW_MLR, XUP_MLR);
				
                hQLvsQR[mod][lay][fi] = new TH2D(QLvsQR.c_str(), title_QLvsQR.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC, NBINS_QDC, XLOW_QDC, XUP_QDC);
                hTLvsTR[mod][lay][fi] = new TH2D(TLvsTR.c_str(), title_TLvsTR.c_str(), NBINS_T, XLOW_T, XUP_T, NBINS_T, XLOW_T, XUP_T);                
                hQn1NormvsQn[mod][lay][fi][0] = new TH2D(Qn1NormvsQn_L.c_str(), title_Qn1NormvsQn_L.c_str(), 50, 0, 1, 50, XLOW_QDC, XUP_QDC);
                hQn1NormvsQn[mod][lay][fi][1] = new TH2D(Qn1NormvsQn_R.c_str(), title_Qn1NormvsQn_R.c_str(), 50, 0, 1, 50, XLOW_QDC, XUP_QDC);
                hQn_1NormvsQn[mod][lay][fi][0] = new TH2D(Qn_1NormvsQn_L.c_str(), title_Qn_1NormvsQn_L.c_str(), 50, 0, 1, 50, XLOW_QDC, XUP_QDC);
                hQn_1NormvsQn[mod][lay][fi][1] = new TH2D(Qn_1NormvsQn_R.c_str(), title_Qn_1NormvsQn_R.c_str(), 50, 0, 1, 50, XLOW_QDC, XUP_QDC);
                
                activeAddresses.push_back({mod, lay, fi});
            }     
            
            if(pRaw->getQDCL() > QDC_MIN &&
            sigL->GetAmplitude() < ADC_MAX &&
            sigL->GetTOT() > TOT_MIN &&
            pRaw->getTimeL() > T_MIN)
            {
                auxiliaryCounts[mod][lay][fi][0] = 1;
                auxiliaryQvalues[mod][lay][fi][0] = pRaw->getQDCL();
            }
            
            if(pRaw->getQDCR() > QDC_MIN &&
            sigR->GetAmplitude() < ADC_MAX &&
            sigR->GetTOT() > TOT_MIN &&
            pRaw->getTimeR() > T_MIN)
            {
                auxiliaryCounts[mod][lay][fi][1] = 1;
                auxiliaryQvalues[mod][lay][fi][1] = pRaw->getQDCR();
            }

            if(pRaw->getQDCL() > QDC_MIN &&
            pRaw->getQDCR() > QDC_MIN &&
            sigL->GetAmplitude() < ADC_MAX &&
            sigR->GetAmplitude() < ADC_MAX &&
            sigL->GetTOT() > TOT_MIN &&
            sigR->GetTOT() > TOT_MIN &&
            pRaw->getTimeL() > T_MIN &&
            pRaw->getTimeR() > T_MIN &&
            fabs(pRaw->getTimeL() - pRaw->getTimeR()) < TDIFF_MAX /*&&
            sigL->GetVeto() == 0 &&
            sigR->GetVeto() == 0*/)
            {
                hA[mod][lay][fi][0]->Fill(sigL->GetAmplitude());
                hA[mod][lay][fi][1]->Fill(sigR->GetAmplitude());
                hQ[mod][lay][fi][0]->Fill(pRaw->getQDCL());
                hQ[mod][lay][fi][1]->Fill(pRaw->getQDCR());
                hT[mod][lay][fi][0]->Fill(pRaw->getTimeL());
                hT[mod][lay][fi][1]->Fill(pRaw->getTimeR());
                
                hQAve[mod][lay][fi]->Fill(sqrt(pRaw->getQDCL()*pRaw->getQDCR()));
                hQSum[mod][lay][fi]->Fill(pRaw->getQDCL() + pRaw->getQDCR());
                hMLR[mod][lay][fi]->Fill(log(sqrt(pRaw->getQDCR()/pRaw->getQDCL())));
                
                hQLvsQR[mod][lay][fi]->Fill(pRaw->getQDCL(), pRaw->getQDCR());
                hTLvsTR[mod][lay][fi]->Fill(pRaw->getTimeL(), pRaw->getTimeR());
                
                ++mult;
                
            }	

            if(pRaw->getQDCL() > QDC_MIN &&
            pRaw->getQDCR() > QDC_MIN &&
            sigL->GetAmplitude() < ADC_MAX &&
            sigR->GetAmplitude() < ADC_MAX &&
            sigL->GetTOT() > TOT_MIN &&
            sigR->GetTOT() > TOT_MIN &&
            pRaw->getTimeL() > T_MIN &&
            pRaw->getTimeR() > T_MIN /*&&
            sigL->GetVeto() == 0 &&
            sigR->GetVeto() == 0*/)
            {
                hTDif[mod][lay][fi]->Fill(pRaw->getTimeL() - pRaw->getTimeR() );
            }
			
		}
// 		cout << auxiliaryCounts[mod][lay][0][0] << " " << auxiliaryCounts[mod][lay][1][0] << " " << auxiliaryCounts[mod][lay][2][0] << " " << auxiliaryCounts[mod][lay][3][0] << " " << auxiliaryCounts[mod][lay][4][0] << " " << auxiliaryCounts[mod][lay][5][0] << " " << auxiliaryCounts[mod][lay][6][0] << " " << auxiliaryCounts[mod][lay][7][0] << endl;
        for(int fibno = FIRST_ACTIVE_FIBER; fibno <= LAST_ACTIVE_FIBER; ++fibno){
            if(auxiliaryCounts[FIXED_MODULE][FIXED_LAYER][fibno][0] == 1 && auxiliaryCounts[FIXED_MODULE][FIXED_LAYER][fibno][1] == 1){
                nCoincidences++;
                if(auxiliaryQvalues[FIXED_MODULE][FIXED_LAYER][fibno][0]+auxiliaryQvalues[FIXED_MODULE][FIXED_LAYER][fibno][1] > maxDepositionFiber) maxDepositionFiber = fibno;
            }
        }

        if(nCoincidences >= 1){
//             for (int fibno = FIRST_ACTIVE_FIBER; fibno <= LAST_ACTIVE_FIBER; ++fibno){
//                 
//             }
            for(int fibside = 0; fibside <=1; fibside++){
                for (int fibno = FIRST_ACTIVE_FIBER; fibno <= LAST_ACTIVE_FIBER; ++fibno){
                    if(auxiliaryCounts[FIXED_MODULE][FIXED_LAYER][fibno][fibside] != 0){
                        if(fibno==maxDepositionFiber){
                            hCrosstalk[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][fibside]->SetBinContent(fibno+1,hCrosstalk[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][fibside]->GetBinContent(fibno+1)+1);
                        }
                        if(fibno < maxDepositionFiber){
                            for(int k=fibno;k<maxDepositionFiber;k++){
                                if(auxiliaryCounts[FIXED_MODULE][FIXED_LAYER][k][fibside] == 0) gapBetweenNeighbours = 1;
                            }
                            if(gapBetweenNeighbours == 0){                                                                              hCrosstalk[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][fibside]->SetBinContent(fibno+1,hCrosstalk[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][fibside]->GetBinContent(fibno+1)+1);
                            if(maxDepositionFiber > 0)
                                hQn_1NormvsQn[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][fibside]->Fill(auxiliaryQvalues[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber+1][fibside]/auxiliaryQvalues[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][fibside], auxiliaryQvalues[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][fibside]);
                            }
                        }
                        gapBetweenNeighbours = 0;
                        if(fibno > maxDepositionFiber){
                            for(int k=maxDepositionFiber;k<fibno;k++){
                                if(auxiliaryCounts[FIXED_MODULE][FIXED_LAYER][k][fibside] == 0) gapBetweenNeighbours = 1;
                            }
                            if(gapBetweenNeighbours == 0){
                                hCrosstalk[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][fibside]->SetBinContent(fibno+1,hCrosstalk[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][fibside]->GetBinContent(fibno+1)+1); //does this method skips a small fraction of entries in comparison to Fill()?
                                if(maxDepositionFiber < 7)
                                    hQn1NormvsQn[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][fibside]->Fill(auxiliaryQvalues[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber+1][fibside]/auxiliaryQvalues[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][fibside], auxiliaryQvalues[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][fibside]);
                            }
                        }
                        gapBetweenNeighbours = 0;
                    

    //----- Use this instead of the above if you don't want to use cut for no gap between neighbours

    //                 if(auxiliaryCounts[FIXED_MODULE][FIXED_LAYER][fibno][0] != 0){
    //                     hCrosstalk[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][0]->SetBinContent(fibno+1,hCrosstalk[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][0]->GetBinContent(fibno+1)+1);
    //                     
    //                 if(auxiliaryCounts[FIXED_MODULE][FIXED_LAYER][fibno][1] != 0){
    //                     hCrosstalk[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][1]->SetBinContent(fibno+1,hCrosstalk[FIXED_MODULE][FIXED_LAYER][maxDepositionFiber][1]->GetBinContent(fibno+1)+1);

                    }
                }
            }
        }
        
        if(nCoincidences > N_ACTIVE_FIBERS || nCoincidences < 0) cout << "wrong coincidence number" << endl;
        
		hFiberMult->Fill(mult);
        nCoincidences = 0;
		mult = 0;
        //auxiliaryCounts[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES] = {};
        maxDepositionFiber = -1;
		loop->nextEvent();
        memset(auxiliaryCounts, 0, sizeof(auxiliaryCounts));
        memset(auxiliaryQvalues, 0, sizeof(auxiliaryQvalues));
	}
	
	std::cout << "\n\nLoop entries: " << nLoop << std::endl;

	path.ReplaceAll(".root", "_HISTOS.root");
	TFile *output = new TFile(path,"RECREATE");
	output->cd();

    for(int i = FIRST_ACTIVE_FIBER; i <=LAST_ACTIVE_FIBER; i++){
//         cout << "n entries in -1st neighbour L: " << hCrosstalk[FIXED_MODULE][FIXED_LAYER][i][0]->GetBinContent(hCrosstalk[FIXED_MODULE][FIXED_LAYER][i][0]->GetMaximumBin()-1) << endl;
        cout /*<< "n entries in -1st neighbour R: "*/ << hCrosstalk[FIXED_MODULE][FIXED_LAYER][i][1]->GetBinContent(hCrosstalk[FIXED_MODULE][FIXED_LAYER][i][1]->GetMaximumBin()-1)/hCrosstalk[FIXED_MODULE][FIXED_LAYER][i][1]->GetBinContent(hCrosstalk[FIXED_MODULE][FIXED_LAYER][i][1]->GetMaximumBin()) << endl;
//         cout << "n entries in +1st neighbour L: " << hCrosstalk[FIXED_MODULE][FIXED_LAYER][i][0]->GetBinContent(hCrosstalk[FIXED_MODULE][FIXED_LAYER][i][0]->GetMaximumBin()+1) << endl;
//         cout << "n entries in +1st neighbour R: " << hCrosstalk[FIXED_MODULE][FIXED_LAYER][i][1]->GetBinContent(hCrosstalk[FIXED_MODULE][FIXED_LAYER][i][1]->GetMaximumBin()+1) << endl;
    }
    
     for(int i = FIRST_ACTIVE_FIBER; i <=LAST_ACTIVE_FIBER; i++){
         cout /*<< "n entries in +1st neighbour R: " */<< hCrosstalk[FIXED_MODULE][FIXED_LAYER][i][1]->GetBinContent(hCrosstalk[FIXED_MODULE][FIXED_LAYER][i][1]->GetMaximumBin()+1)/hCrosstalk[FIXED_MODULE][FIXED_LAYER][i][1]->GetBinContent(hCrosstalk[FIXED_MODULE][FIXED_LAYER][i][1]->GetMaximumBin()) << endl;
    }
    
//----- canvas definition    
    TCanvas * canQDC = new TCanvas("QDC", "QDC", XCANVAS, YCANVAS); 
    canQDC->SetLeftMargin(0.2);
	canQDC->DivideSquare(activeAddresses.size());

//     TCanvas * canT = new TCanvas("T", "T", XCANVAS, YCANVAS); 
//     canT->SetLeftMargin(0.2); 
//     canT->DivideSquare(activeAddresses.size());
    
    TCanvas * canCrosstalk = new TCanvas("canCrosstalk", "canCrosstalk", XCANVAS, YCANVAS); 
    canCrosstalk->SetLeftMargin(0.2); 
    canCrosstalk->DivideSquare(activeAddresses.size());
    
    TCanvas * canQn1NormvsQn_L = new TCanvas("canQn1NormvsQn_L", "canQn1NormvsQn_L", XCANVAS, YCANVAS); 
    canQn1NormvsQn_L->SetLeftMargin(0.2); 
    canQn1NormvsQn_L->DivideSquare(activeAddresses.size());
    
    TCanvas * canQn1NormvsQn_R = new TCanvas("canQn1NormvsQn_R", "canQn1NormvsQn_R", XCANVAS, YCANVAS); 
    canQn1NormvsQn_R->SetLeftMargin(0.2); 
    canQn1NormvsQn_R->DivideSquare(activeAddresses.size());
    
    TCanvas * canQn_1NormvsQn_L = new TCanvas("canQn_1NormvsQn_L", "canQn_1NormvsQn_L", XCANVAS, YCANVAS); 
    canQn_1NormvsQn_L->SetLeftMargin(0.2); 
    canQn_1NormvsQn_L->DivideSquare(activeAddresses.size());
    
    TCanvas * canQn_1NormvsQn_R = new TCanvas("canQn_1NormvsQn_R", "canQn_1NormvsQn_R", XCANVAS, YCANVAS); 
    canQn_1NormvsQn_R->SetLeftMargin(0.2); 
    canQn_1NormvsQn_R->DivideSquare(activeAddresses.size());

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

//----- canvas plotting  
    for(int i=1; i< activeAddresses.size(); i++){ //changed from i=0 bc fiber 0 invalid in this series
        m=activeAddresses[i].iMod;
        l=activeAddresses[i].iLay;
        f=activeAddresses[i].iFib;

//         if(!output->GetListOfKeys()->Contains(Form("Module%d", m))) output->mkdir(Form("Module%d/", m));
//         output->cd(Form("Module%d/", m));
//         if(!gDirectory->GetListOfKeys()->Contains(Form("Layer%d", l))) output->mkdir(Form("Module%d/Layer%d", m, l));
//         output->cd(Form("Module%d/Layer%d", m, l));
//         if(!gDirectory->GetListOfKeys()->Contains(Form("Fiber%d", f))) output->mkdir(Form("Module%d/Layer%d/Fiber%d",m, l, f));
//         output->cd(Form("Module%d/Layer%d/Fiber%d",m, l, f));   
//         hQ[m][l][f][0]->SetDirectory(gDirectory);
//         hQ[m][l][f][0]->Write();
// 
//         hQ[m][l][f][0]->SetDirectory(gDirectory);
//         hQ[m][l][f][0]->Write();
//         hQ[m][l][f][1]->SetDirectory(gDirectory);
//         hQ[m][l][f][1]->Write();
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
        
        canQDC->cd(i+1);
        gPad->SetGrid(1, 1);
        hQ[m][l][f][0]->SetLineColor(colL);
        hQ[m][l][f][1]->SetLineColor(colR);

        hQ[m][l][f][0]->Draw();
        hQ[m][l][f][1]->Draw("same");
        PaintBin(hQ[m][l][f][0], f+1, kBlue);
        PaintBin(hQ[m][l][f][1], f+1, kBlue);
        hQ[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hQ[m][l][f][0]->GetBinContent(hQ[m][l][f][0]->GetMaximumBin()), hQ[m][l][f][1]->GetBinContent(hQ[m][l][f][1]->GetMaximumBin())) + 50);
        format_h_1D(hQ[m][l][f][0]);
        format_h_1D(hQ[m][l][f][1]);
        
//         canT->cd(i+1);
//         gPad->SetGrid(1, 1);
//         hT[m][l][f][0]->SetLineColor(colL);
//         hT[m][l][f][1]->SetLineColor(colR);
//         hT[m][l][f][0]->Draw();
//         hT[m][l][f][1]->Draw("same");
//         hT[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hT[m][l][f][0]->GetBinContent(hT[m][l][f][0]->GetMaximumBin()), hT[m][l][f][1]->GetBinContent(hT[m][l][f][1]->GetMaximumBin())) + 50);
//         format_h_1D(hT[m][l][f][0]);
//         format_h_1D(hT[m][l][f][1]);
        
        canCrosstalk->cd(i); //changed from i+1
        gPad->SetGrid(1, 1);
        gPad->SetLogy();
        hCrosstalk[m][l][f][0]->SetLineColor(colL);
        hCrosstalk[m][l][f][1]->SetLineColor(colR);
        hCrosstalk[m][l][f][0]->Scale(1./hCrosstalk[m][l][f][0]->GetBinContent(hCrosstalk[m][l][f][0]->GetMaximumBin()));
        hCrosstalk[m][l][f][1]->Scale(1./hCrosstalk[m][l][f][1]->GetBinContent(hCrosstalk[m][l][f][1]->GetMaximumBin()));
        hCrosstalk[m][l][f][0]->GetYaxis()->SetRangeUser(1e-5*2, 2);
        hCrosstalk[m][l][f][0]->Draw();
        hCrosstalk[m][l][f][1]->Draw("same");
        PaintBin(hCrosstalk[m][l][f][0], f+1, kBlue);
//         PaintBin(hCrosstalk[m][l][f][1], f+1, kBlue);
//         hCrosstalk[m][l][f][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hCrosstalk[m][l][f][0]->GetBinContent(hCrosstalk[m][l][f][0]->GetMaximumBin()), hCrosstalk[m][l][f][1]->GetBinContent(hCrosstalk[m][l][f][1]->GetMaximumBin())) + 50);
        format_h_1D(hCrosstalk[m][l][f][0]);
        format_h_1D(hCrosstalk[m][l][f][1]); 

        
        canQn1NormvsQn_L->cd(i); //changed from i+1
        gPad->SetGrid(1, 1);
        hQn1NormvsQn[m][l][f][0]->Draw("colz");
        format_h_1D(hQn1NormvsQn[m][l][f][0]);
        
        canQn1NormvsQn_R->cd(i); //changed from i+1
        gPad->SetGrid(1, 1);
        hQn1NormvsQn[m][l][f][1]->Draw("colz");
        format_h_1D(hQn1NormvsQn[m][l][f][1]);
        
        canQn_1NormvsQn_L->cd(i); //changed from i+1
        gPad->SetGrid(1, 1);
        hQn_1NormvsQn[m][l][f][0]->Draw("colz");
        format_h_1D(hQn_1NormvsQn[m][l][f][0]);
            
        canQn_1NormvsQn_R->cd(i); //changed from i+1
        gPad->SetGrid(1, 1);
        hQn_1NormvsQn[m][l][f][1]->Draw("colz");
        format_h_1D(hQn_1NormvsQn[m][l][f][1]);
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
               
    } 
canCrosstalk->cd(8);
        
    auto legend = new TLegend(0.1,0.3,0.9,0.7); // option "C" allows to center the header
    legend->AddEntry(hCrosstalk[FIXED_MODULE][FIXED_LAYER][6][0],"LEFT","l");
    legend->AddEntry(hCrosstalk[FIXED_MODULE][FIXED_LAYER][6][1],"RIGHT","l");
    legend->AddEntry((TObject*)0, "marked bin: coinc. & highest dep. E ", "");
//     legend->AddEntry(hCrosstalk[FIXED_MODULE][FIXED_LAYER][6][0],"Coupled to fiber","l");
//     legend->AddEntry(hCrosstalk[m][l][f][0],"1st neighbour","p");
//     legend->AddEntry(gEventsLN2,"2nd neighbour","p");
    legend->Draw();
    
    
	canMult->cd();  
    gPad->SetGrid(1, 1);
	hFiberMult->Draw();
 	format_h_1D(hFiberMult);
    hFiberMult->GetXaxis()->SetRange(1, activeAddresses.size()+1);
	hFiberMult->GetXaxis()->SetTitle("# fibers");
	hFiberMult->GetYaxis()->SetTitle("Counts");
    
// 	output->mkdir("Canvases");
// 	output->cd("Canvases");
//     
// 	canQDC->Write();
//     canT->Write();
// 	canTDif->Write();
//     canQAve->Write();
// 	canQSum->Write();
//     canMLR->Write();
// 	canQLvsQR->Write(); 
// 	canTLvsTR->Write();    
//     canMult->Write();
    
//    output->Close();

//----- writing canvas to file    
//     TString pdfname = ConstructPdfName(path);
//     canQDC->Print(Form("%s(",pdfname.Data()), "pdf");
//     canT->Print(pdfname,"pdf");
//     canTDif->Print(pdfname,"pdf");
//     canQAve->Print(pdfname,"pdf");
//     canQSum->Print(pdfname,"pdf");
//     canMLR->Print(pdfname,"pdf");
//     canQLvsQR->Print(pdfname,"pdf");
//     canTLvsTR->Print(pdfname,"pdf");
//     canMult->Print(Form("%s)",pdfname.Data()),"pdf");

    auto end = std::chrono::system_clock::now();   
    std::chrono::duration<double> time = end-start;
    std::cout << "Time: " << time.count() << " s\n";
    
	return 0;
}
