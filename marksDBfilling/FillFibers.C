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
    sqlite3_stmt *statement;
    sqlite3_prepare_v2(database, query, -1, &statement, nullptr);
    sqlite3_step(statement);
    sqlite3_finalize(statement);
}
void FillFibers(UInt_t seriesIDvalue, const char * database_name = "prototype_db.sqlite3") {
    gSystem->Load("libsqlite3.so");  
    sqlite3 *database;
    int status = sqlite3_open(database_name, &database);
    if(status!=0){
        std::cerr << "##### Could not access data base!" << std::endl;
        return 0;
    }
    const char * query = Form("select measurement_id, Measurements.series_id, n_modules, n_layers_per_module,  n_fibers_per_layer from Measurements join Series where Measurements.series_id=Series.series_id and Measurements.series_id=%d", seriesIDvalue); //changed from series_id>17;
    std::vector<std::map<std::string, std::string> > results = Query(database, query);
    for(int i=0; i < results.size(); ++i) {
        printf("%d / %zu\n", i, results.size() );
        for(UInt_t moduleID=0; moduleID < std::stoi(results[i]["n_modules"]); ++moduleID){ //this loop is added
            for(UInt_t layer=0; layer < std::stoi(results[i]["n_layers_per_module"]); ++layer){
                for(UInt_t fiber=0; fiber < std::stoi(results[i]["n_fibers_per_layer"]); ++fiber){
                    //                 Exec(database, Form("INSERT INTO Fibers (measurement_id, series_id, module, layer, fiber) VALUES (%s, %s, %s, %d, %d)", results[i]["measurement_id"].c_str(), results[i]["series_id"].c_str(), results[i]["n_modules"].c_str(), layer, fiber) );
                    Exec(database, Form("INSERT INTO Fibers (measurement_id, series_id, module, layer, fiber) VALUES (%s, %s, %d, %d, %d)", results[i]["measurement_id"].c_str(), results[i]["series_id"].c_str(), moduleID, layer, fiber) );
                }
            }                  
        }                            
    }
}
