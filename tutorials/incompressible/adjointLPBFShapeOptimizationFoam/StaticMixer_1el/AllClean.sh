# !/bin/bash

echo -e "   - removing old results"
rm -rf 1* 2* 3* 4* 5* 6* 7* 8* 9*
rm -r processor*

echo -e "   - Removing folder 0"
rm -r 0

echo -e "   - Removing old mesh"
rm -r constant/polyMesh

echo -e "   - Removing dynamicCode"
rm -r dynamicCode

echo -e "   - Removing old log files"
rm log.*