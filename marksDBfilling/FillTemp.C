#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sqlite3.h>
#include <vector>

#include "TSystem.h"
#include "TString.h"

bool FillTemp(int seriesNo){

std::cout << "test" << std::endl;
    
gSystem->Load("libsqlite3.so");    
    
TString query;
sqlite3 *database;
sqlite3_stmt * statement;
int status;

TString tempfile;
std::vector <TString> names;
std::vector <int> start;
std::vector <int> stop;

TString sensorID;
TString celc;
int time;
double temperature;
    
//----- opening database 

// TString databasename = "/home/kasia/data/DB/ScintFib_2.db";
TString databasename = "/scratch/gccb/kasia/data/DB/ScintFib_2.db";
status = sqlite3_open(databasename, &database);

if(status!=0){
  std::cerr << "Couldn't open data base!" << std::endl;
  return false;
}

//----- reading series details
query = Form("SELECT TEMP_FILE FROM SERIES WHERE SERIES_ID = %i", seriesNo);
status = sqlite3_prepare_v2(database, query, -1, &statement, nullptr);

if(status!=0){
  std::cerr << "Problem executing query #1" << std::endl;
  return false;
}

while((status = sqlite3_step(statement)) == SQLITE_ROW){
  const unsigned char *tempname_char = sqlite3_column_text(statement, 0);
  tempfile = std::string(reinterpret_cast<const char*>(tempname_char));
}

status = sqlite3_finalize(statement);
    
std::cout << "Temperature log for series " << seriesNo << " is: ";
std::cout << tempfile << std::endl;

//----- reading measurement details
query = Form("SELECT MEASUREMENT_NAME, START_TIME, STOP_TIME FROM MEASUREMENT WHERE SERIES_ID = %i", seriesNo);    
status = sqlite3_prepare_v2(database, query, -1, &statement, nullptr);

if(status!=0){
  std::cerr << "Problem executing query #2" << std::endl;
  return false;
}

while((status = sqlite3_step(statement))==SQLITE_ROW){
  const unsigned char *name_char = sqlite3_column_text(statement, 0);
  names.push_back(std::string(reinterpret_cast <const char*> (name_char)));
  start.push_back(sqlite3_column_int(statement, 1));
  stop.push_back(sqlite3_column_int(statement, 2));
  std::cout << "\tReading data base..." << std::endl;
}

status = sqlite3_finalize(statement);

int n = names.size();
std::cout << "List of measurements in series ";
std::cout << seriesNo << ": " << std::endl;

for(int i=0; i<n; i++){
  std::cout << names[i] << "\t" << start[i];
  std::cout << "\t" << stop[i] << std::endl;
}

//----- filling temperatures

TString tempfile_full = "/scratch/gccb/kasia/data/temp_logs/"+tempfile;
//ifstream input(tempfile_full);
ifstream input;

for(int i=0; i<n; i++){
  //input.seekg(0);
  //std::streampos pos = input.tellg();
  //std::cout << "Position in file: " << pos << std::endl;
  input.open(tempfile_full);
  if(!input.is_open()){
    std::cerr << "Couldn't open temperature log!" << std::endl;
    return false;
  }
  while(input.good()){
    input >> sensorID >> time >> temperature >> celc;
      if(time>start[i] && time<stop[i]){
        query = Form("INSERT INTO TEMPERATURES (SERIES_ID, SENSOR_ID, MEASUREMENT_NAME, TIME, TEMPERATURE) VALUES (%i, '%s', '%s', %lld, %.2f);", seriesNo, sensorID.Data(), names[i].Data(), (long long) time, temperature);
        status = sqlite3_prepare_v2(database, query, -1, &statement, nullptr);
        //std::cout << query << std::endl;
        if(status!=0){
          std::cerr << "Problem executing query #3" << std::endl;
          return false;
        }
        status = sqlite3_step(statement);
      }
   }
   std::cout << "Filling temperatures for file ";
   std::cout << names[i] << " done!" << std::endl;
   input.close();
}

sqlite3_finalize(statement);

//input.close();

return true;
}
