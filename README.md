# phd-thesis-sasse

The primary purpose of this repository is to serve as supplementary material to the Dissertation "Adjoint Topology Optimisation of Polymer Melt Flow Channels Producible by Additive Manufacturing" by Jana Sasse.

![Proposed changes for a static mixer that is optimised for minimal pressure loss](/DemoStaticMixer.png)

## Contents
+ Adjoint topology optimisation solver (based on the ``adjointShapeOptimizationFoam`` solver by *Othmer et al.*)
  + Consideration of laser powder bed fusion (LPBF) manufacturing constraints
    + Assumes positive x-direction as print direction
    + Mesh needs to be cartesian between layers but can be non-cartesian within a layer
  + Optimisation objectives for
    + Minimal pressure drop
    + Target temperature at the outlet
    + Target concentration at the outlet
+ Carreau-WLF viscosity model (optimised for polymer melt flows)
+ Tutorial for optimisation of an additively manufacturable 1D static mixer

Please note that this solver and the tutorial were only tested for [OpenFOAM-6](https://openfoam.org/release/6/).

## Citation
Please cite the following papers in any publication for which you find this solver useful.
+ Sasse, J.: *Adjoint Topology Optimisation of Polymer Melt Flow Channels Producible by Additive Manufacturing*. RWTH Aachen University, Dissertation (in preparation)
+ Sasse, J.; Schön, M.; Hopmann, C.: Static Mixers Producible by Additive Manufacturing: Novel Rapid Automatic Optimisation and Practical Evaluation. *Polymers* 14 (2022) 21, p. 4646
+ Sasse, J.; Hopmann, C.: Static mixers producible by additive manufacturing: Operating point specific geometries through automatic optimisation. *Proceedings of the 38th International Conference of the Polymer Processing Society (PPS‐38)*. St. Gallen, CH, 2024
+ Sasse, J.; Hopmann, C.: Automatic optimisation of additively manufactured static mixers under consideration of design restrictions. *Proceedings of the 32nd International Colloquium on Plastics Technology*. Aachen, DE, 2024

## Contact
+ Contributor: Jana Sasse
+ Email:       jana.sasse.aachen@gmail.com

## License
Distributed using the GNU General Public License (GPL), version 3; see the LICENSE file for details.
