#include <iostream>

#include "SCategoryManager.h"
#include "SLoop.h"
#include "SFibersRawCluster.h"

#define N_LAYERS_PER_MODULE 7
#define N_FIBERS_PER_LAYER 55

#define N_HITMAP_TYPES 2
#define N_THRESHOLDS 1

#define QDC_MIN 0
// #define QDC_MAX 7000
// #define QDC_MAX 587
#define QDC_MAX 834

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
    
    TString hitMapName[N_HITMAP_TYPES][N_THRESHOLDS] = {nullptr}; 
    TString enDepoName[N_HITMAP_TYPES][N_THRESHOLDS] = {nullptr}; 
    TH2D * hitMap[N_HITMAP_TYPES][N_THRESHOLDS] = {nullptr};
    TH2D * enDepo[N_HITMAP_TYPES][N_THRESHOLDS] = {nullptr};
    std::vector<char*> type = {"AE", "NS"};
//     std::vector<int> threshold = {0,500,1000,1500}; //EXPERIMENT [keV]
//     std::vector<int> threshold = {0,710,1362,2013}; //SIMULATION [Photon Counts]; corresponds to 0,500,1000,1500 keV
//     std::vector<int> threshold = {451}; //just for the cut around 511
    std::vector<int> threshold = {640}; //just for the cut around 511 - SIMULATION
    for(int i=0; i<N_HITMAP_TYPES; i++){
        for(int j=0; j<N_THRESHOLDS; j++){
            hitMapName[i][j] = Form("hFiberHitMap_%s_t%04d",type[i],threshold[j]);
            enDepoName[i][j] = Form("hEnergyDepoMap_%s_t%04d",type[i],threshold[j]);
            hitMap[i][j] = new TH2D(hitMapName[i][j], hitMapName[i][j], N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
            enDepo[i][j] = new TH2D(enDepoName[i][j], enDepoName[i][j], N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
//             std::cout << hitMapName[i][j] << std::endl;
//             std::cout << enDepoName[i][j] << std::endl;
        }
    }
// printf("Q_L: %s\n", Q_L.Data());
    
    
    
    std::vector<std::string> str_path;
    for(int i=0; i < path.size(); i++){
        str_path.push_back(std::string(path[i]));
    }
    
	SLoop * loop = new SLoop();
	loop->addFiles(str_path);
	loop->setInput({});
    
    //if needed, change the category below to any other category valid in sifi-framework (also add appropriate header and change class object names accordingly)
	SCategory * pCatRawClus = SCategoryManager::getCategory(SCategory::CatFibersRawClus);
    
    TH2D * hFiberHitMap_AE_t0_ERR = new TH2D("hFiberHitMap_AE_t0_ERR", "hFiberHitMap_AE_t0_ERR", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);

    

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
            
            for(int j=0; j<N_THRESHOLDS; j++)
            {
                if(pRawClus->getQDCL() > threshold[j] && pRawClus->getQDCR() > threshold[j] && pRawClus->getQDCL() < QDC_MAX && pRawClus->getQDCR() < QDC_MAX) //added upper threshold
//                 if(pRawClus->getQDCL() > QDC_MIN && pRawClus->getQDCR() > QDC_MIN && pRawClus->getQDCL() < QDC_MAX && pRawClus->getQDCR() < QDC_MAX) //added upper threshold
                {
                    qav=sqrt(pRawClus->getQDCL()*pRawClus->getQDCR());
                    hitMap[0][j]->Fill(fi, lay);
                    enDepo[0][j]->Fill(fi, lay, qav);
                    if(pRawClus->getFiberClusterLabel()>0 && pRawClus->getFiberClusterLabel()<3)
                    {
                        hitMap[1][j]->Fill(fi, lay);
                        enDepo[1][j]->Fill(fi, lay, qav);
                    }
                }
            }
        }
    }

    int binCount=0;
    for(int i=0; i<N_LAYERS_PER_MODULE; i++)
    {
        for(int j=0; j<N_FIBERS_PER_LAYER; j++)
        {
            globalBin=hitMap[0][0]->GetBin(j+1, i+1);
//             std::cout << "globalBin " << globalBin << std::endl;
            binCount = hitMap[0][0]->GetBinContent(globalBin);
//             std::cout << "binCount " << binCount << std::endl;
            if(binCount!=0)
            hFiberHitMap_AE_t0_ERR->SetBinContent(globalBin,1/sqrt(binCount));
            else hFiberHitMap_AE_t0_ERR->SetBinContent(globalBin,0);
        }
    }
    
    int xChannel, yChannel;
    std::ifstream istream;
    istream.open("/scratch/gccb/magda/DeadFibers.txt", std::ios::in);
    if (!istream.is_open())
    {
        std::cerr << "Error! Could not open input file!" << std::endl;
        return false;
    }
//     std::string csvLine;
    istream.seekg(std::ios_base::beg); //set pointer to beginning
    while(!istream.eof()){
        istream >> yChannel >> xChannel;
    //     std::cout << xChannel << " " <<  yChannel << std::endl;
        globalBin=hitMap[0][0]->GetBin(xChannel+1, yChannel+1);

        hFiberHitMap_AE_t0_ERR->SetBinContent(globalBin, 0);
        
        for(int i=0; i<N_HITMAP_TYPES; i++){
            for(int j=0; j<N_THRESHOLDS; j++){
                hitMap[i][j]->SetBinContent(globalBin, 0);
                enDepo[i][j]->SetBinContent(globalBin, 0);
            }
        }
    }
    


	std::cout << "\n\nLoop entries: " << nLoop << std::endl;

	path[ii].ReplaceAll(".root", "_HITMAPS_SIM.root");
	TFile *output = new TFile(path[ii],"RECREATE");
	output->cd();
    
    for(int i=0; i<N_HITMAP_TYPES; i++){
        for(int j=0; j<N_THRESHOLDS; j++){
            hitMap[i][j]->Write();
            enDepo[i][j]->Write();
        }
    }

    
    TCanvas * c1 = new TCanvas("c1", "c1", 1600, 800);
    
    c1->cd();
    hFiberHitMap_AE_t0_ERR->Draw("colz"); 




    auto end = std::chrono::system_clock::now();   
    std::chrono::duration<double> time = end-start;
    std::cout << "Time: " << time.count() << " s\n";
    
	return 0;
}
