#include <iostream>
#include <fstream>
#include <stdlib.h>

void citiroc_spectroscopy_reader(){
    std::ifstream istream;
    istream.open("/scratch1/gccb/data/202207/citiroc/Run662_list.txt", std::ios::in);
//     istream.open("/home/magda/DATA_CITIROC/Run316_list.txt", std::ios::in);
    if (!istream.is_open())
    {
        std::cerr << "Error! Could not open input file!" << std::endl;
        return false;
    }

    std::string csvLine;
    std::string tmp;
    std::string str_tstamp_us, str_trgID, str_brd, str_ch, str_q_lg, str_q_hg;
    float tstamp_us = 0.;
    float q_hg = 0.;
    int trgID = 0;
    int brd = 0; 
    int ch = 0;
    std::string q_lg;
    int not_spaces=0;

    istream.seekg(std::ios_base::beg);
    for (int i=0; i<9;i++) getline(istream, tmp); //skip header
    
    while(!istream.eof()){
        getline(istream, csvLine);
        for(int j=0; j<csvLine.length(); j++){
            if(csvLine[j]!=' ') not_spaces++;
        }
        if(not_spaces <= 10){
            std::stringstream stream(csvLine);
            stream >> str_brd >> str_ch >> str_q_lg >> str_q_hg;
            brd = std::stoi(str_brd);
            ch = std::stoi(str_ch);
            q_hg = std::stof(str_q_hg);
            q_lg = str_q_lg;
        cout <<  brd << " " << ch << " " << q_lg << " " << q_hg << endl;
        }
        if(not_spaces > 10){
            std::stringstream stream(csvLine);
            stream >> str_tstamp_us >> str_trgID >> str_brd >> str_ch >> str_q_lg >> str_q_hg;
            tstamp_us = std::stof(str_tstamp_us);
            trgID = std::stoi(str_trgID);
            brd = std::stoi(str_brd);
            ch = std::stoi(str_ch);
            q_hg = std::stof(str_q_hg);
            q_lg = str_q_lg;
        cout << tstamp_us << " " <<  trgID <<" " <<  brd << " " << ch << " " << q_lg << " " << q_hg << endl;
        }
        not_spaces=0;
    }
}
