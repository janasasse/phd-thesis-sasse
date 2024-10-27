# Tutorial adjointLPBFShapeOptimizationFoam

The primary purpose of this tutorial is to demonstrate the usage of the solver developed alongside the dissertation "Adjoint Topology Optimisation of Polymer Melt Flow Channels Producible by Additive Manufacturing" by Jana Sasse.

## Contents

Please note that this tutorial was only tested for [OpenFOAM-6](https://openfoam.org/release/6/).

+ Initial geometry (files in ``STL`` folder)
  + Assumes build in direction of extrusion (positive x-direction)
  + Additively manufacturable 1D static mixer (``initialGeometry.stl``)
  + Includes scaffold structure (``scaffoldStructure.stl``) and outer cylinder to support application-specific manufacturing constraints
+ Optimisation objectives
  + Optimisation for minimal pressure loss for 500 steps in ``./AllRun.sh``
  + Optimisation for maximal thermal or material homogeneity are available for right choice of weights
+ Operating point
  + Polymer material is a PE-HD modelled with Carreau-WLF
  + Mixer material has thermal properties of steel
  + Throughput corresponds to roughly 20 kg/h
  + Inlet thermal inhomogeneity corresponds to hot core profile
  + Inlet material inhomogeneity corresponds to half profile

![Optimisation for minimal pressure loss](/tutorials/incompressible/adjointLPBFShapeOptimizationFoam/StaticMixer_1el/optPtLoss.png)

## How to run

+ Use ``./AllClean.sh`` to delete files created in previous run
+ Type ``./AllRun.sh`` in your bash terminal
  + Initialises case
  + Computes solution for geometry before optimisation (parallel run)
  + Performs optimisation for minimal pressure loss for 500 steps (serial run)
+ Recommendation: calculate a 'clean' solution of state variables with optimised geometry for evaluation (especially after optimising for thermal or material mixing)
  + copy ``0.orig`` to ``0``
  + copy ``alpha`` from optimised result folder to folder ``0``
  + in ``constant/optimizationProperties``:
    + set ``initialization`` flag to ``yes``
    + set ``wP``, ``wT``, ``wC`` to ``0``
  + adjust system/controlDict such that solver calculates solution until convergence
  + run application ``adjointLPBFShapeOptimizationFoam`` (can be run in parallel)

## How to view and interpret results

View results in [Paraview](https://www.paraview.org/), e.g., using the ``StaticMixer_1el.foam`` file provided

The optimised geometry is embedded in scalarField ``alpha`` and can by visualised by one of the following filters
+ Threshold of ``alpha`` between ``1e+09`` and ``2e+09`` (will give most exact result)
+ Clip (scalar) of ``alpha`` at ``1e+09`` (make sure upper half of that threshold is visualised; will give smoothest visualisation)

![Recommended visualisation techniques](/tutorials/incompressible/adjointLPBFShapeOptimizationFoam/StaticMixer_1el/visualisationTechniques.png)

``alphaDiff`` can be used to highlight the proposed direction of change at the interface. In this case, red highlights areas where algorithm wants to add material to the static mixer, while blue highlights areas where algorithm wants to remove material from the mixer.

![Highlighting the proposed changes using alphaDiff](/tutorials/incompressible/adjointLPBFShapeOptimizationFoam/StaticMixer_1el/alphaDiff.png)

The optimised geometry can be exported with the following steps
1. Clip (scalar) of ``alpha`` at ``1e+09`` (make sure upper half of that threshold is visualised)
2. Extract surface
3. Save data as STL
