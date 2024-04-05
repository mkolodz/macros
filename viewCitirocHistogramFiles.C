#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLine.h>
#include <TString.h>
#include <TTree.h>

#define file_len 8192
Bool_t viewCitirocHistogramFiles(TString path)
{
    if(!path.Contains("/") || !path.BeginsWith("/"))
    {
        std::cout << "Error! The functions needs the filename including the full absolute path..." << std::endl;
        std::abort();
    }

    std::ifstream istream;
    istream.open(path, std::ios::in);
    if (!istream.is_open()) //this is important because otherwise there is no warning if file does not exist!
    {
        std::cerr << "Error in SKSSource::open()! Could not open input file!" << std::endl;
        return false;
    }
    std::string csvLine;
    istream.seekg(std::ios_base::beg); //set pointer to beginning

    int energy_channel[file_len]={0};
    int bincount[file_len]={0};
    int line_counter=0;

    while(line_counter!=file_len){

        istream >> csvLine; //std::stoi(getline(istream, csvLine)); 
        bincount[line_counter] = std::stoi(csvLine);
//         cout << energy_channel << " " << bincount << endl; //prints first line
        energy_channel[line_counter]=line_counter;
        line_counter++;
        
    }
    
    TGraph *g= new TGraph(file_len, energy_channel, bincount);
    TCanvas *c1 = new TCanvas("c1","Graph Draw Options",200,10,600,400);
    g->Draw("AP");
    
//     c1->Close();
    
    
    return 0;
}
