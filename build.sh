g++ -I/scratch/gccb/software/include -L/scratch/gccb/software/lib -L/scratch/gccb/software/lib64 open_fitparams_for_full_series.C -lScintillatingFibers -lSiFi -lFibers -lCore $(root-config --cflags --libs) -o open_fitparams.o
#g++ -I/scratch/gccb/software/include -L/scratch/gccb/software/lib -L/scratch/gccb/software/lib64 copy_fitparams.C -lScintillatingFibers -lSiFi -lFibers -lCore $(root-config --cflags --libs) -o copy_fitparams.o

