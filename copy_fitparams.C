#include <iostream>
#include <fstream>
#include <string> 
#include <stdlib.h>
#include <cstring>  

int main(int argc, char **argv)
{
  if(argc!=2)
  {
    std::cout << "to run type: ./copy_fitparams.o path/to/logfile.txt" << std::endl;
    return 0;
  }
  
  char* log_name = argv[1];
  
  int nletters = strlen(log_name);
  
  int istop = -1;
  
  for (int i = nletters; i > 0; i--)
  {
      if (log_name[i] == '/')
      {
          istop = i;
          break;
      }
  }

  int istart = 0;
   
  if (istop == -1)
  {
      std::cerr << "##### Error in main()!" << std::endl;
      std::cerr << "Cannot interpret log name!" << std::endl;
      return 0;
  }

  std::string data_path = std::string(&log_name[istart], &log_name[istop + 1]);
  std::cout << "data path: " << data_path << std::endl;
  
  std::ifstream log(log_name); 
  
  if(!log.is_open())
  {
   std::cerr << "Couldn't open log file!" << std::endl;  
   std::cerr << log_name << std::endl;
   return 0;
  }
  
  //std::string log_name_folder = std::string(&log_name[istart], &log_name[istop]);
  std::string dummy;
  std::string meas_name;
  std::string meas_path;
  int status = 0;
  std::string log_name_folder = std::string(&log_name[nletters-20], &log_name[nletters-4]);
  std::cout << log_name_folder << std::endl;
  
  
  while(log.good())
  {
   getline(log,dummy);
   getline(log,dummy);
   getline(log,dummy);
   log >> meas_name;
   meas_path = std::string(data_path) + "../" + meas_name + "/";
   getline(log,dummy);
   getline(log,dummy);
   getline(log,dummy);
   
   char command_3[900] = "cp ";
   strcat(command_3, meas_path.c_str());
   strcat(command_3, "fitparams.out /scratch/gccb/data/fitparams_copy/");
   //strcat(command_3, meas_path.c_str());
   strcat(command_3, log_name_folder.c_str());
   strcat(command_3, "/");
   strcat(command_3, meas_name.c_str());
   strcat(command_3, "_fitparams.out");

   std::cout << command_3 << std::endl;
   puts(command_3);
   status = system(command_3);
   std::cout << "Status: " << status << std::endl;
  }

  return 1;
}
