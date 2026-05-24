#!/bin/bash
echo Starting Computation
nohup matlab -nodesktop -nosplash -nodisplay < RunOptimisation.m > logfile0.txt 2> err0.txt &
nohup matlab -nodesktop -nosplash -nodisplay < RunSlave.m > logfile1.txt 2> err1.txt &
nohup matlab -nodesktop -nosplash -nodisplay < RunSlave.m > logfile2.txt 2> err2.txt &
nohup matlab -nodesktop -nosplash -nodisplay < RunSlave.m > logfile3.txt 2> err3.txt &
