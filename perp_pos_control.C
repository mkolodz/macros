#include "SCategoryManager.h"
#include "SLoop.h"
#include "SDDSamples.h"
#include "SFibersRaw.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH2I.h"
#include "THStack.h"
#include "TTree.h"
#include "TString.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TMath.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <RtypesCore.h>

#define NBINS 100
#define N_MODULES 1
#define N_LAYERS 4
#define N_FIBERS 16
#define N_SIDES 2
#define XLOW_ADC 0
#define XUP_ADC 3000
#define XLOW_QDC -500
#define XUP_QDC 5000
#define XLOW_Qavg 0
#define XUP_Qavg 3000
#define XLOW_Time 0
#define XUP_Time 1000
#define XCANVAS 1600
#define YCANVAS 1000
#define ARR_LEN 50000
#define XUP_QDC_REL 1
#define XLOW_QDC_REL 0
using namespace std;

void format_h_31(TH1 * h) {
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.0);
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.2);
    h->GetYaxis()->SetMaxDigits(3);
    h->GetZaxis()->SetLabelSize(0.04);
    gPad->SetLeftMargin(0.2);
    gPad->SetTopMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.23);
    gPad->SetFrameLineWidth(1); 
    gPad->SetTicky(1);
    gPad->SetGrid(1,1);
    gStyle->SetTitleFontSize(0.07);
}

int perp_pos_control(std::string name="sifi_results.root")
// int rawTimes(TString name="/scratch/gccb/data/202109/2021_09_10_14_48/sifi_results.root")
//int rawTimes(std::vector<Double_t>& entriesL, std::string name="/scratch/gccb/data/202109/2021_09_10_14_48/sifi_results.root")
{
    //std::vector<Double_t>& entriesL,
    //gStyle->SetOptStat(1111);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(57);
    //gStyle->SetLabelSize(2);
    
    SLoop * loop = new SLoop();
    loop->addFile(name);
    loop->setInput({});
    
    TFile *output = new TFile("perp_pos_control.root","RECREATE");
    SCategory * pCatRaw = SCategoryManager::getCategory(SCategory::CatFibersRaw);
    SCategory * pCatSample = SCategoryManager::getCategory(SCategory::CatDDSamples);
    Int_t mod, lay, fi, side;
    Int_t n_active = 0;
    TString histNames[512];

    TH1D * hA[N_MODULES][N_LAYERS][N_FIBERS][N_SIDES];
    TH1D * hQ[N_MODULES][N_LAYERS][N_FIBERS][N_SIDES];
    TH1D * hT[N_MODULES][N_LAYERS][N_FIBERS][N_SIDES];
    TH1D * hTavg[N_MODULES][N_LAYERS][N_FIBERS];
    TH1D * hQavg[N_MODULES][N_LAYERS][N_FIBERS];
    TH1D * hQsum[N_MODULES][N_LAYERS][N_FIBERS];
    TH2D * hQlQr[N_MODULES][N_LAYERS][N_FIBERS];
    TH2D * hTlTr[N_MODULES][N_LAYERS][N_FIBERS];
    TH1D * hTDiff[N_MODULES][N_LAYERS][N_FIBERS];
    

    int n = loop->getEntries();
    for (int i = 0; i < n; ++i)
    {
        size_t nn = pCatRaw->getEntries();
        for (uint j = 0; j < nn; ++j)
        {
            SDDSamples* samples = (SDDSamples*)pCatSample->getObject(j);
            SFibersRaw * pRaw = (SFibersRaw *)pCatRaw->getObject(j);
            SDDSignal* sigL = (SDDSignal*)samples->getSignalL();
            SDDSignal* sigR = (SDDSignal*)samples->getSignalR();
            pRaw->getAddress(mod, lay, fi);
            //std::cout << "HERE 1" << std::endl;
            TString A_L = Form("A_M%dL%dF%dS0",mod,lay,fi);
            TString A_R = Form("A_M%dL%dF%dS1",mod,lay,fi);
            TString Q_L = Form("Q_M%01dL%01dF%02dS0",mod,lay,fi);
            TString Q_R = Form("Q_M%01dL%01dF%02dS1",mod,lay,fi);
            TString T_L = Form("T_M%01dL%01dF%02dS0",mod,lay,fi);
            TString T_R = Form("T_M%01dL%01dF%02dS1",mod,lay,fi);
            TString Qavg = Form("Qavg_M%01dL%01dF%02d",mod,lay,fi);
            TString Tavg = Form("Tavg_M%01dL%01dF%02d",mod,lay,fi);
            TString Qsum = Form("Qsum_M%01dL%01dF%02d",mod,lay,fi);
            TString QlQr = Form("QlQr_M%01dL%01dF%02d",mod,lay,fi);
            TString TlTr = Form("TlTr_M%01dL%01dF%02d",mod,lay,fi);
            TString Tdiff = Form("Tdiff_M%01dL%01dF%02d",mod,lay,fi);
                   //     std::cout << "HERE 2" << std::endl;

            if(!gDirectory->FindObject(Q_L)){
                hA[mod][lay][fi][0] = new TH1D(A_L, Form("%s; A_{l} [mV]; Counts", A_L.Data()), NBINS, 0, 100);
                hA[mod][lay][fi][1] = new TH1D(A_R, Form("%s; A_{r} [mV]; Counts", A_R.Data()), NBINS, 0, 100);
                hQ[mod][lay][fi][0] = new TH1D(Q_L, Form("%s; Q_{l} [nVs]; Counts", Q_L.Data()), NBINS, XLOW_QDC, XUP_QDC);
                hQ[mod][lay][fi][1] = new TH1D(Q_R, Form("%s; Q_{r} [nVs]; Counts", Q_R.Data()), NBINS, XLOW_QDC, XUP_QDC);
                hT[mod][lay][fi][0] = new TH1D(T_L, Form("%s; T_{l} [ns]; Counts", T_L.Data()), NBINS, -200,1000);
                hT[mod][lay][fi][1] = new TH1D(T_R, Form("%s; T_{r} [ns]; Counts", T_R.Data()), NBINS, -200,1000);
                hQavg[mod][lay][fi] = new TH1D(Qavg, Form("%s; Q_{avg} [nVs]; Counts", Qavg.Data()), NBINS, XLOW_QDC, XUP_QDC);
                hTavg[mod][lay][fi] = new TH1D(Tavg, Form("%s; T_{avg} [ns]; Counts", Tavg.Data()), NBINS, XLOW_ADC, 1000);
                hQsum[mod][lay][fi] = new TH1D(Qsum, Form("%s; Q_{l} + Q_{r} [nVs]; Counts", Qsum.Data()), NBINS, XLOW_ADC, XUP_ADC);
                hQlQr[mod][lay][fi] = new TH2D(QlQr, Form("%s; Q_{l} [nVs]; Q_{r} [nVs]", QlQr.Data()), NBINS, XLOW_QDC, XUP_QDC, NBINS, XLOW_QDC, XUP_QDC);
                hTlTr[mod][lay][fi] = new TH2D(TlTr, Form("%s; T_{l} [ns]; T_{r} [ns]", TlTr.Data()), NBINS, XLOW_Time, XUP_Time, NBINS, XLOW_Time, XUP_Time);
                hTDiff[mod][lay][fi] = new TH1D(Tdiff, Form("%s; T_{diff} [ns]; Counts", Tavg.Data()), NBINS, -100, 100);
                
                histNames[n_active] = Q_L;
                n_active++;
            }
                      //  std::cout << "HERE 3" << std::endl;

            hA[mod][lay][fi][0]->Fill(sigL->GetAmplitude());
            hA[mod][lay][fi][1]->Fill(sigR->GetAmplitude());
                     //   std::cout << "HERE 4" << std::endl;

            if(pRaw->getQDCL()!=-100){
                hQ[mod][lay][fi][0]->Fill(pRaw->getQDCL());
            }
            if(pRaw->getQDCR()!=-100){
                hQ[mod][lay][fi][1]->Fill(pRaw->getQDCR());
            }
            if(pRaw->getTimeL()!=-100){
                hT[mod][lay][fi][0]->Fill(pRaw->getTimeL());
            }
            if(pRaw->getTimeR()!=-100){
                hT[mod][lay][fi][1]->Fill(pRaw->getTimeR());
            }
            if(pRaw->getQDCR()!=-100 && pRaw->getQDCL()!=-100){
                hQavg[mod][lay][fi]->Fill(sqrt(pRaw->getQDCL()*pRaw->getQDCR()));
                hQsum[mod][lay][fi]->Fill(pRaw->getQDCL()+pRaw->getQDCR());
                hQlQr[mod][lay][fi]->Fill(pRaw->getQDCL(), pRaw->getQDCR());
            }
            if(pRaw->getTimeL()!=-100 && pRaw->getTimeR()!=-100){
                hTavg[mod][lay][fi]->Fill(sqrt(pRaw->getTimeL()*pRaw->getTimeR()));
                hTDiff[mod][lay][fi]->Fill(pRaw->getTimeL()-pRaw->getTimeR());
                hTlTr[mod][lay][fi]->Fill(pRaw->getTimeL(),pRaw->getTimeR());
            }
        }
        loop->nextEvent();
    }
    std::cout << "loop_entries: " << n << std::endl;

   // output->cd();
//     for(Int_t m = 0; m < N_MODULES; m++)
//     {
//         output->mkdir(Form("Module%d/", m));
//         output->cd(Form("Module%d/", m));
//         for(Int_t l = 0; l < N_LAYERS; l++)
//         {
//             output->mkdir(Form("Module%d/Layer%d", m, l));
//             for(Int_t f = 0; f < N_FIBERS; f++)
//             {
//                 output->mkdir(Form("Module%d/Layer%d/Fiber%d", m,l,f));
//                 output->cd(Form("Module%d/Layer%d/Fiber%d", m,l,f));
//                 hQ[m][l][f][0]->SetDirectory(gDirectory);
//                 hQ[m][l][f][0]->Write();
//                 hQ[m][l][f][1]->SetDirectory(gDirectory);
//                 hQ[m][l][f][1]->Write();
//                 hA[m][l][f][0]->SetDirectory(gDirectory);
//                 hA[m][l][f][0]->Write();
//                 hA[m][l][f][1]->SetDirectory(gDirectory);
//                 hA[m][l][f][1]->Write();
//                 hT[m][l][f][0]->SetDirectory(gDirectory);
//                 hT[m][l][f][0]->Write();
//                 hT[m][l][f][1]->SetDirectory(gDirectory);
//                 hT[m][l][f][1]->Write();
//                 hQavg[m][l][f]->SetDirectory(gDirectory);
//                 hQavg[m][l][f]->Write();
//                 hQsum[m][l][f]->SetDirectory(gDirectory);
//                 hQsum[m][l][f]->Write();
//                 hTavg[m][l][f]->SetDirectory(gDirectory);
//                 hTavg[m][l][f]->Write();
//                 hQlQr[m][l][f]->SetDirectory(gDirectory);
//                 hQlQr[m][l][f]->Write();
//                 hTDiff[m][l][f]->SetDirectory(gDirectory);
//                 hTDiff[m][l][f]->Write();
//                 hTlTr[m][l][f]->SetDirectory(gDirectory);
//                 hTlTr[m][l][f]->Write();
//             }
//         }
//     }    

    TCanvas * canQDCL = new TCanvas("QDCL", "QDCL", XCANVAS, YCANVAS); 
    canQDCL->DivideSquare(n_active);
    
    TCanvas * canQDCR = new TCanvas("canQDCR", "QDCR", XCANVAS, YCANVAS); 
    canQDCR->DivideSquare(n_active);
    
    TCanvas * canADCL = new TCanvas("ADCL", "ADCL", XCANVAS, YCANVAS); 
    canADCL->DivideSquare(n_active);
    
    TCanvas * canADCR = new TCanvas("ADCR", "ADCR", XCANVAS, YCANVAS); 
    canADCR->DivideSquare(n_active);
    
    TCanvas * canTimeL = new TCanvas("TimeL", "TimeL", XCANVAS, YCANVAS); 
    canTimeL->DivideSquare(n_active);
    
    TCanvas * canTimeR = new TCanvas("TimeR", "TimeR", XCANVAS, YCANVAS); 
    canTimeR->DivideSquare(n_active);
    
    TCanvas * canQavg = new TCanvas("Qavg", "Qavg", XCANVAS, YCANVAS); 
    canQavg->DivideSquare(n_active);
    
    TCanvas * canTavg = new TCanvas("Tavg", "Tavg", XCANVAS, YCANVAS); 
    canTavg->DivideSquare(n_active);
    
    TCanvas * canQsum = new TCanvas("Qsum", "Qsum", XCANVAS, YCANVAS); 
    canQsum->DivideSquare(n_active);
    
    TCanvas * canQlQr = new TCanvas("QlQr", "QlQr", XCANVAS, YCANVAS); 
    canQlQr->DivideSquare(n_active);
    
    TCanvas * canTlTr = new TCanvas("TlTr", "TlTr", XCANVAS, YCANVAS); 
    canTlTr->DivideSquare(n_active);
    
    TCanvas * canMult = new TCanvas("Mult", "Mult",600,400); 
    
    TCanvas * canTime = new TCanvas("Time", "Time", XCANVAS, YCANVAS); 
    
    canTime->Divide(3,1);
    canTime->cd(1);  
    canTime->SetLeftMargin(0.2);
    hT[0][0][0][0]->Draw();
    format_h_31(hT[0][0][0][0]);
    canTime->cd(2);
    hT[0][0][0][1]->Draw();
    format_h_31(hT[0][0][0][1]);
    canTime->cd(3);
    hTDiff[0][0][0]->Draw();
    format_h_31(hTDiff[0][0][0]);
    
//      for(Int_t i = 0; i < N_MODULES; i++)
//      {
//         for(Int_t j = 0; j < N_LAYERS; j++)
//         {
//             for(Int_t k = 0; k < N_FIBERS; k++)
//             {
                canQDCL->cd(1);  
                canQDCL->SetLeftMargin(0.2);
                hQ[0][j][k][0]->Draw();
                format_h_31(hQ[i][j][k][0]);
                canQDCL->cd(2);  
                canQDCL->SetLeftMargin(0.2);
                hQ[i][j][k][0]->Draw();
                format_h_31(hQ[i][j][k][0]);
                canQDCL->cd(3);  
                canQDCL->SetLeftMargin(0.2);
                hQ[i][j][k][0]->Draw();
                format_h_31(hQ[i][j][k][0]);
                canQDCL->cd(4);  
                canQDCL->SetLeftMargin(0.2);
                hQ[i][j][k][0]->Draw();
                format_h_31(hQ[i][j][k][0]);
                
                canQDCR->cd(i*N_LAYERS*N_FIBERS+j*N_FIBERS+k+1);  
                canQDCR->SetLeftMargin(0.2);
                hQ[i][j][k][1]->Draw();
                format_h_31(hQ[i][j][k][1]);
/*                
                canADCL->cd(i*N_LAYERS*N_FIBERS+j*N_FIBERS+k+1);  
                canADCL->SetLeftMargin(0.2);
                hA[i][j][k][0]->Draw();
                format_h_31(hA[i][j][k][0]);
                
                canADCR->cd(i*N_LAYERS*N_FIBERS+j*N_FIBERS+k+1);  
                canADCR->SetLeftMargin(0.2);
                hA[i][j][k][1]->Draw();
                format_h_31(hA[i][j][k][1]);
                
                canTimeL->cd(i*N_LAYERS*N_FIBERS+j*N_FIBERS+k+1);  
                canTimeL->SetLeftMargin(0.2);
                hT[i][j][k][0]->Draw();
                format_h_31(hT[i][j][k][0]);
                
                canTimeR->cd(i*N_LAYERS*N_FIBERS+j*N_FIBERS+k+1);  
                canTimeR->SetLeftMargin(0.2);
                hT[i][j][k][1]->Draw();
                format_h_31(hT[i][j][k][1]);
                
                canQavg->cd(i*N_LAYERS*N_FIBERS+j*N_FIBERS+k+1);  
                canQavg->SetLeftMargin(0.2);
                hQavg[i][j][k]->Draw();
                format_h_31(hQavg[i][j][k]);
                
                canTavg->cd(i*N_LAYERS*N_FIBERS+j*N_FIBERS+k+1);  
                canTavg->SetLeftMargin(0.2);
                hTavg[i][j][k]->Draw();
                format_h_31(hTavg[i][j][k]);

                canQsum->cd(i*N_LAYERS*N_FIBERS+j*N_FIBERS+k+1);  
                canQsum->SetLeftMargin(0.2);
                hQsum[i][j][k]->Draw();
                format_h_31(hQsum[i][j][k]);
                
                canQlQr->cd(i*N_LAYERS*N_FIBERS+j*N_FIBERS+k+1);  
                canQlQr->SetLeftMargin(0.2);
                gStyle->SetPalette(kBird);
                hQlQr[i][j][k]->Draw("COLZ");
                format_h_31(hQlQr[i][j][k]);
                
                canTlTr->cd(i*N_LAYERS*N_FIBERS+j*N_FIBERS+k+1);  
                canTlTr->SetLeftMargin(0.2);
                gStyle->SetPalette(kBird);
                hTlTr[i][j][k]->Draw("COLZ");
                format_h_31(hTlTr[i][j][k]);*/

            }
        }
    }

//     Int_t counter=0;
//     for (Int_t i=0;i<sizeof(histNames)/sizeof(histNames[0]); i++){
//         if(histNames[i]!="") 
//         {
//             counter++;
//         }
//     }
// 
//    std::cout <<  "# of histograms per side: " << counter << std::endl;
//     
//     canMult->cd();  
//     hFiberMult->Draw();
//     format_h_31(hFiberMult);
//     hFiberMult->GetXaxis()->SetTitle("N_fibers");
//     hFiberMult->GetYaxis()->SetTitle("Counts");
// 
//     output->mkdir("Canvases");
//     output->cd("Canvases");
//     canQDCL->Write();
//     canQDCR->Write();
//     canADCL->Write();    
//     canADCR->Write();
//     canTimeL->Write();
//     canTimeR->Write();
//     canTavg->Write();
//     canQavg->Write();
//     canQsum->Write();
//     canQlQr->Write();   
//     canTlTr->Write();
   
     return 0;
}


// int main() {
//     rawTimes();
// }
