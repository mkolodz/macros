#include <iostream>

#include "SCategoryManager.h"
#include "SLoop.h"
#include "SFibersRawCluster.h"

//----- prototype geometry 4to1
#define N_LAYERS_PER_MODULE 7
#define N_FIBERS_PER_LAYER 55

//----- cuts
#define QDC_MIN 0

//how to  run: 
//root -l
//.x rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel.C({"/path/to/sifi_tree1.root","/path/to/sifi_tree2.root"})

int rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel(std::vector<TString> path)
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
    
    //if needed, change the category below to any other category valid in sifi-framework (also add appropriate header and change class object names accordingly)
	SCategory * pCatRawClus = SCategoryManager::getCategory(SCategory::CatFibersRawClus);
    
    TH2D * hFiberHitMap_AE_t0 = new TH2D("hFiberHitMap_AE_t0", "hFiberHitMap_AE_t0", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hFiberHitMap_AE_t500 = new TH2D("hFiberHitMap_AE_t500", "hFiberHitMap_AE_t500", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hFiberHitMap_AE_t1000 = new TH2D("hFiberHitMap_AE_t1000", "hFiberHitMap_AE_t1000", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hFiberHitMap_AE_t1500 = new TH2D("hFiberHitMap_AE_t1500", "hFiberHitMap_AE_t1500", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    
    TH2D * hFiberHitMap_NS_t0 = new TH2D("hFiberHitMap_NS_t0", "hFiberHitMap_NS_t0", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hFiberHitMap_NS_t500 = new TH2D("hFiberHitMap_NS_t500", "hFiberHitMap_NS_t500", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hFiberHitMap_NS_t1000 = new TH2D("hFiberHitMap_NS_t1000", "hFiberHitMap_NS_t1000", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hFiberHitMap_NS_t1500 = new TH2D("hFiberHitMap_NS_t1500", "hFiberHitMap_NS_t1500", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    
    TH2D * hEnergyDepoMap_AE_t0 = new TH2D("hEnergyDepoMap_AE_t0", "hEnergyDepoMap_AE_t0", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hEnergyDepoMap_AE_t500 = new TH2D("hEnergyDepoMap_AE_t500", "hEnergyDepoMap_AE_t500", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hEnergyDepoMap_AE_t1000 = new TH2D("hEnergyDepoMap_AE_t1000", "hEnergyDepoMap_AE_t1000", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hEnergyDepoMap_AE_t1500 = new TH2D("hEnergyDepoMap_AE_t1500", "hEnergyDepoMap_AE_t1500", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    
    TH2D * hEnergyDepoMap_NS_t0 = new TH2D("hEnergyDepoMap_NS_t0", "hEnergyDepoMap_NS_t0", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hEnergyDepoMap_NS_t500 = new TH2D("hEnergyDepoMap_NS_t500", "hEnergyDepoMap_NS_t500", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hEnergyDepoMap_NS_t1000 = new TH2D("hEnergyDepoMap_NS_t1000", "hEnergyDepoMap_NS_t1000", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hEnergyDepoMap_NS_t1500 = new TH2D("hEnergyDepoMap_NS_t1500", "hEnergyDepoMap_NS_t1500", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);

    int QDC_UP=7000;
	Int_t mod, lay, fi;
    int globalBin=0;
    double qav=0;
    SFibersRawCluster * pRawClus;
 	int nLoop = loop->getEntries();
    
    for(int i=0; i<nLoop; i++)
    {
        loop->getEvent(i);
        size_t nCat = pCatRawClus->getEntries();
		for (uint j = 0; j < nCat; j++)
        {
            pRawClus = (SFibersRawCluster *)pCatRawClus->getObject(j);
            pRawClus->getAddress(mod, lay, fi);
            if(pRawClus->getQDCL() > QDC_MIN && pRawClus->getQDCR() > QDC_MIN)
            {
                qav=sqrt(pRawClus->getQDCL()*pRawClus->getQDCR());
                hFiberHitMap_AE_t0->Fill(fi, lay);
                hEnergyDepoMap_AE_t0->Fill(fi, lay, qav); //geometric mean of QDCL and QDCR
                
                if(pRawClus->getFiberClusterLabel()>0 && pRawClus->getFiberClusterLabel()<3)
                {
                    hFiberHitMap_NS_t0->Fill(fi, lay);
                    hEnergyDepoMap_NS_t0->Fill(fi, lay, qav);
                }
            }
            
            if(pRawClus->getQDCL() > 500 && pRawClus->getQDCR() > 500)
            { 
                qav=sqrt(pRawClus->getQDCL()*pRawClus->getQDCR());
                hFiberHitMap_AE_t500->Fill(fi, lay);
                hEnergyDepoMap_AE_t500->Fill(fi, lay, qav); //geometric mean of QDCL and QDCR
                
                if(pRawClus->getFiberClusterLabel()>0 && pRawClus->getFiberClusterLabel()<3)
                {
                    hFiberHitMap_NS_t500->Fill(fi, lay);
                    hEnergyDepoMap_NS_t500->Fill(fi, lay, qav);
                }
            }
            if(pRawClus->getQDCL() > 1000 && pRawClus->getQDCR() > 1000)
            { 
                qav=sqrt(pRawClus->getQDCL()*pRawClus->getQDCR());
                hFiberHitMap_AE_t1000->Fill(fi, lay);
                hEnergyDepoMap_AE_t1000->Fill(fi, lay, qav); //geometric mean of QDCL and QDCR
                
                if(pRawClus->getFiberClusterLabel()>0 && pRawClus->getFiberClusterLabel()<3)
                {
                    hFiberHitMap_NS_t1000->Fill(fi, lay);
                    hEnergyDepoMap_NS_t1000->Fill(fi, lay, qav);
                }
            }
                
            if(pRawClus->getQDCL() > 1500 && pRawClus->getQDCR() > 1500)
            { 
                qav=sqrt(pRawClus->getQDCL()*pRawClus->getQDCR());
                hFiberHitMap_AE_t1500->Fill(fi, lay);
                hEnergyDepoMap_AE_t1500->Fill(fi, lay, qav); //geometric mean of QDCL and QDCR
                
                if(pRawClus->getFiberClusterLabel()>0 && pRawClus->getFiberClusterLabel()<3)
                {
                    hFiberHitMap_NS_t1500->Fill(fi, lay);
                    hEnergyDepoMap_NS_t1500->Fill(fi, lay, qav);
                }
            }
        }
    }

    int xChannel, yChannel;
    std::ifstream istream;
    istream.open("DeadFibers.txt", std::ios::in);
    if (!istream.is_open())
    {
        std::cerr << "Error! Could not open input file!" << std::endl;
        return false;
    }
//     std::string csvLine;
    istream.seekg(std::ios_base::beg); //set pointer to beginning
    while(!istream.eof()){
    istream >> yChannel >> xChannel;
    std::cout << xChannel << " " <<  yChannel << std::endl;
    globalBin=hFiberHitMap_AE_t0->GetBin(xChannel+1, yChannel+1);

    hFiberHitMap_AE_t0->SetBinContent(globalBin, 0);
    hFiberHitMap_AE_t500->SetBinContent(globalBin, 0);
    hFiberHitMap_AE_t1000->SetBinContent(globalBin, 0);
    hFiberHitMap_AE_t1500->SetBinContent(globalBin, 0);

    hFiberHitMap_NS_t0->SetBinContent(globalBin, 0);
    hFiberHitMap_NS_t500->SetBinContent(globalBin, 0);
    hFiberHitMap_NS_t1000->SetBinContent(globalBin, 0);
    hFiberHitMap_NS_t1500->SetBinContent(globalBin, 0);

    hEnergyDepoMap_AE_t0->SetBinContent(globalBin, 0);
    hEnergyDepoMap_AE_t500->SetBinContent(globalBin, 0);
    hEnergyDepoMap_AE_t1000->SetBinContent(globalBin, 0);
    hEnergyDepoMap_AE_t1500->SetBinContent(globalBin, 0);

    hEnergyDepoMap_NS_t0->SetBinContent(globalBin, 0);
    hEnergyDepoMap_NS_t500->SetBinContent(globalBin, 0);
    hEnergyDepoMap_NS_t1000->SetBinContent(globalBin, 0);
    hEnergyDepoMap_NS_t1500->SetBinContent(globalBin, 0);
    }

	std::cout << "\n\nLoop entries: " << nLoop << std::endl;

	path[ii].ReplaceAll(".root", "_HITMAPS_chooseFiberLabel.root");
	TFile *output = new TFile(path[ii],"RECREATE");
	output->cd();
    hFiberHitMap_AE_t0->Write();
    hFiberHitMap_AE_t500->Write();
    hFiberHitMap_AE_t1000->Write();
    hFiberHitMap_AE_t1500->Write();

    hFiberHitMap_NS_t0->Write();
    hFiberHitMap_NS_t500->Write();
    hFiberHitMap_NS_t1000->Write();
    hFiberHitMap_NS_t1500->Write();
                                     
    hEnergyDepoMap_AE_t0->Write();
    hEnergyDepoMap_AE_t500->Write();
    hEnergyDepoMap_AE_t1000->Write();
    hEnergyDepoMap_AE_t1500->Write();
                                 
    hEnergyDepoMap_NS_t0->Write();
    hEnergyDepoMap_NS_t500->Write();
    hEnergyDepoMap_NS_t1000->Write();
    hEnergyDepoMap_NS_t1500->Write();

    
//     TCanvas * c1 = new TCanvas("c1", "c1", 1600, 800);
//     c1->Divide(4,4);
//     
//     c1->cd(1);
// hFiberHitMap_AE_t0->Draw("colz"); 
// c1->cd(2);
// hFiberHitMap_AE_t500->Draw("colz"); 
// c1->cd(3);
// hFiberHitMap_AE_t1000->Draw("colz");
// c1->cd(4);
// hFiberHitMap_AE_t1500->Draw("colz");
// 
// c1->cd(5);
// hFiberHitMap_NS_t0->Draw("colz"); 
// c1->cd(6);
// hFiberHitMap_NS_t500->Draw("colz"); 
// c1->cd(7);
// hFiberHitMap_NS_t1000->Draw("colz");
// c1->cd(8);
// hFiberHitMap_NS_t1500->Draw("colz");
// 
// 
// c1->cd(9);
// hEnergyDepoMap_AE_t0->Draw("colz");
// c1->cd(10);
// hEnergyDepoMap_AE_t500->Draw("colz");
// c1->cd(11);
// hEnergyDepoMap_AE_t1000->Draw("colz"); 
// c1->cd(12);
// hEnergyDepoMap_AE_t1500->Draw("colz"); 
// 
// c1->cd(13);
// hEnergyDepoMap_NS_t0->Draw("colz");
// c1->cd(14);
// hEnergyDepoMap_NS_t500->Draw("colz");
// c1->cd(15);
// hEnergyDepoMap_NS_t1000->Draw("colz"); 
// c1->cd(16);
// hEnergyDepoMap_NS_t1500->Draw("colz"); 


    auto end = std::chrono::system_clock::now();   
    std::chrono::duration<double> time = end-start;
    std::cout << "Time: " << time.count() << " s\n";
    
	return 0;
}
