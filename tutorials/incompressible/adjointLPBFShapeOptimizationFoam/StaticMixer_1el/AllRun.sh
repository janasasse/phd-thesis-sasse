#!/bin/bash

echo -e "\n   This is an algorithm that performs adjoint-based topology optimisation"
echo -e " under consideration of the design restrictions from additive manufacturing\n"

# Number of procs for parallel computation
nProcs=8

# Number of steps performed during optimisation
optimisationSteps=500

# Definition of optimization objective
wP=1
wT=0
wC=0

foamDictionary constant/optimizationProperties -entry wP -set "wP [0 0 0 0 0 0 0] $wP"
foamDictionary constant/optimizationProperties -entry wT -set "wT [0 0 0 0 0 0 0] $wT"
foamDictionary constant/optimizationProperties -entry wC -set "wC [0 0 0 0 0 0 0] $wC"

# Meshing and IB initialization
echo -e "\n   - Set up case\n"
  echo -e "     ++ Copy 0.orig to 0"
    rm -r 0
    cp -r 0.orig 0
  echo -e "     ++ Create initial mesh"
    blockMesh &> log.setupCase
  echo -e "     ++ Rotate mesh for build direction / direction of extrusion in +x"
    transformPoints -rotate '((1 0 0)(0 0 -1))' &>> log.setupCase
  echo -e "     ++ Initialize geometry with setFields\n"
    setFields &>> log.setupCase

# Initialization run
echo -e "\n   - Compute initial primal and adjoint CFD solution (in parallel with $nProcs procs)\n"
  echo -e "     ++ Prepare files"
    foamDictionary constant/optimizationProperties -entry initialization -set yes
    foamDictionary system/decomposeParDict -entry numberOfSubdomains -set $nProcs
    foamDictionary system/controlDict -entry startFrom -set latestTime
    foamDictionary system/controlDict -entry startTime -set 0
    foamDictionary system/controlDict -entry endTime -set 10000
    foamDictionary system/controlDict -entry writeInterval -set 10000
  echo -e "     ++ Decompose domain"
    decomposePar &> log.initializationRun
    wait
  echo -e "     ++ Compute SIMPLE solution"
    mpirun -np $nProcs adjointLPBFShapeOptimizationFoam -parallel &>> log.initializationRun
    wait
  echo -e "     ++ Reconstruct domain and clean up\n"
    reconstructPar &>> log.initializationRun
    wait

# Optimization Run
echo -e "\n   - Prepare optimization run\n"
  rm -r processor*
  endInitRun=`grep 'SIMPLE solution converged' log.initializationRun | cut -d' ' -f5`
  endOptRun=$(echo "scale=3; $endInitRun + $optimisationSteps" | bc -l)
  foamDictionary constant/optimizationProperties -entry initialization -set no
  foamDictionary system/controlDict -entry startFrom -set latestTime
  foamDictionary system/controlDict -entry startTime -set $endInitRun
  foamDictionary system/controlDict -entry endTime -set $endOptRun
  foamDictionary system/controlDict -entry writeInterval -set 100
echo -e "\n   - Perform optimization (in serial) for $optimisationSteps steps"
echo -e "     optimization objective: wP = $wP, wT = $wT, wC = $wC\n"
  adjointLPBFShapeOptimizationFoam &> log.optimizationRun

echo -e "\n Finished optimisation!\n"
echo -e " The optimised geometry is contained in the scalarField alpha and can be visualized in Paraview"
echo -e " Check README.md for more details\n"
echo -e " Note that it might be beneficial to recompute a 'clean' CFD solution using the new alpha field before evaluation\n"