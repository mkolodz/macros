#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <ctime>
#include <algorithm>
#include <string>
int nlc_for_calibration(std::string input){
    std::ifstream istream;
    istream.open(input, std::ios::in);
    if (!istream.is_open()) //this is important because otherwise there is no warning if file does not exist!
    {
        std::cerr << "Error: Could not open input file!" << std::endl;
        return false;
    }
    float p0, p1, p2, p3;
    p0=8.0;
    p1=1.04676;
    p2=1.02734;
    p3=0.31909;	
    std::string name;
    double mean, e_mean, sigma, e_sigma;
    double corr_mean, corr_e_mean, corr_sigma, corr_e_sigma;
    std::string csvLine;
    
    ofstream output;
    output.open("rescaled_calib_factors.txt");
    output.clear();
    output << "Writing this to a file.\n";

    
    istream.seekg(std::ios_base::beg);
    getline(istream, csvLine);
    //while (!istream.eof())
    while(true)
    {
//         getline(istream, csvLine);
        istream >> name >> mean >> e_mean >> sigma >> e_sigma;
        if(istream.eof()) break;
        corr_mean = (p0*pow(p1, pow(mean, p2))+p3*mean-p0); 
        corr_e_mean = (p0*pow(p1, pow(e_mean, p2))+p3*e_mean-p0);
        corr_sigma = (p0*pow(p1, pow(sigma, p2))+p3*sigma-p0);
        corr_e_sigma = (p0*pow(p1, pow(e_sigma, p2))+p3*e_sigma-p0);
        cout << name << " "<< mean << " "<< e_mean << " "<< sigma << " "<< e_sigma << " " << endl;
        output << name << " " << corr_mean << " "<< corr_e_mean << " "<< corr_sigma << " "<< corr_e_sigma << " " << endl;

    }
    
    output.close();
                
                
    return 0;
}
