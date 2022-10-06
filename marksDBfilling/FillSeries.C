#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sqlite3.h>
void tokenize(std::string str, std::vector<std::string> &token_v, char delimiter = '/'){
    size_t start = str.find_first_not_of(delimiter), end=start;
    while (start != std::string::npos) {
        // Find next occurence of delimiter
        end = str.find(delimiter, start);
        // Push back the token found into vector
        token_v.push_back(str.substr(start, end-start));
        // Skip all occurences of the delimiter to find new start
        start = str.find_first_not_of(delimiter, end);
    }
}
void findFiberAddress(const char * f, std::vector<std::map<std::string, std::string> > &address) {
    std::ifstream file(f);
    std::string line;
    Bool_t start = kFALSE;
    while(getline(file, line) ) {
        if (line.find("[FibersDDLookupTable]") != std::string::npos) {
            getline(file, line);
            getline(file, line);
            start = kTRUE;
        }
        if(start) {
            std::string::size_type pos;
            while((pos = line.find("  ") ) != std::string::npos) {
                line.replace(pos, 2, " ");
            }
            std::vector<std::string> token;
            tokenize(line, token, ' ');
            std::map<std::string, std::string> mapAddr;
            //consider just the left so we don't repeat
            if(token[5].compare("l") == 0) {
                mapAddr["m"] = token[2];
                mapAddr["l"] = token[3];
                mapAddr["f"] = token[4];
                address.push_back(mapAddr);
            }
        }
    }
}
std::vector<std::map<std::string, std::string> > Query(sqlite3 *database, const char * query) {
    //printf("%s\n", query);
    std::vector<std::map<std::string, std::string> > result;
    sqlite3_stmt *statement;
    sqlite3_prepare_v2(database, query, -1, &statement, nullptr);
    while(sqlite3_step(statement) == SQLITE_ROW) {
        std::map<std::string, std::string> row;
        for(UShort_t i=0; i < sqlite3_data_count(statement); ++i) {
            row[std::string(sqlite3_column_name(statement, i) )] = std::string((char *)sqlite3_column_text(statement, i) );
        }
        result.push_back(row);
    }
    sqlite3_finalize(statement);
    return result;
}
void Exec(sqlite3 *database, const char * query) {
    printf("%s\n", query);
//    sqlite3_stmt *statement;
//    sqlite3_prepare_v2(database, query, -1, &statement, nullptr);
//    sqlite3_step(statement);
//    sqlite3_finalize(statement);
}
void FillSeries(const char * database_name = "prototype_db.sqlite3") {
    gSystem->Load("libsqlite3.so");  
    sqlite3 *database;
    int status = sqlite3_open(database_name, &database);
    if(status!=0){
        std::cerr << "##### Could not access data base!" << std::endl;
        return 0;
    }
    //insert into Series(detector, fiber_mat, fiber_dim, fiber_prod, fiber_wrap, fiber_geometry, n_modules, n_layers_per_module, n_fibers_per_layer, radiosource, photodetector, voltage, coupling_type, coupling_size, daq, configuration, log_file , monitoring_file, pos_unc, mask, elog_id, description) 
    const char * query = "SELECT series_id, log_file FROM Series WHERE series_id IN (87)";
    std::vector<std::map<std::string, std::string> > results = Query(database, query);
    for(int i=0; i < results.size(); ++i) {
        printf("%d / %zu\n", i, results.size() );
        std::vector<std::string> tokenlogs;
        //there are maybe more than one log_file, separated by semicolons
        tokenize(results[i]["log_file"], tokenlogs, ';');
        for(UShort_t k = 0; k < tokenlogs.size(); ++k) {            
            std::ifstream logfile(tokenlogs[k]);
            if(!logfile.is_open() ) {
                std::cout << "Couldn't open log file " << tokenlogs[k] << std::endl;
                return 0;
            }
            std::vector<std::string> tokens;
            //get the directory, in the  YYYYMM format
            tokenize(tokenlogs[k], tokens);
            const char * dir = tokens[3].c_str();

            time_t start_time, stop_time;
            TString meas_name;
            double pos_motor, pos_true;
            std::string dummy;
            while(logfile.good()){
                getline(logfile, dummy);
                logfile >> dummy >> dummy >> start_time;
                start_time += 60;
                logfile >> dummy >> dummy >> stop_time;                
                logfile >> meas_name;
                if(meas_name.IsNull() ) continue;
                logfile >> dummy >> dummy >> pos_motor;
                pos_true = pos_motor + 10;
                getline(logfile, dummy);
                getline(logfile, dummy);
                //get fiber information from FibersDDLookupTable from a Measurement
                std::vector<std::map<std::string, std::string> > vecFiberAddress;
                const char *pathToParams = Form("/scratch1/gccb/data/%s/%s/params.txt", dir, meas_name.Data() );
                findFiberAddress(pathToParams, vecFiberAddress);
                //INSERT Measurements and Fibers into the sqlite3 database
                Exec(database, Form("INSERT INTO Measurements (series_id, datadir, start_time, stop_time, duration, source_pos) VALUES (%s, '/scratch1/gccb/data/%s/%s', %ld, %ld, %ld, '0;%2.f')", results[i]["series_id"].c_str(), dir, meas_name.Data(), start_time, stop_time, stop_time - start_time, pos_true) );
                for(UInt_t j=0; j < vecFiberAddress.size(); ++j) {
                    Exec(database, Form("INSERT INTO Fibers (measurement_id, series_id, module, layer, fiber) VALUES ((SELECT measurement_id FROM Measurements ORDER BY measurement_id DESC LIMIT 1), %s, %s, %s, %s)", results[i]["series_id"].c_str(), vecFiberAddress[j]["m"].c_str(), vecFiberAddress[j]["l"].c_str(), vecFiberAddress[j]["f"].c_str() ) );
                }
            }
        }
    }
}
