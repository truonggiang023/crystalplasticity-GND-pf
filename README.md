# Crystal orientation driven local deformation and damage in dual phase steel with different morphologies and quench-induced residual strains 


This repository contains the source codes developed and implemented in the accompanying manuscript for a computational framework integrating crystal plasticity, geometrically necessary dislocations (GNDs) through thermal expansion FEM, and phase-field damage modeling to investigate the deformation and fracture behavior of dual-phase (DP) steels.

The proposed framework incorporates heterogeneous initial GND distributions induced by martensitic transformation into subsequent mechanical simulations, enabling a quantitative assessment of the coupled effects of crystallographic texture, transformation-induced residual deformation, and martensite morphology on local deformation and damage evolution.

# Code: 
umat2phase.f - the main user-defined material subroutine (UMAT) implementing the constitutive framework for dual-phase (DP) steels in Abaqus. The framework combines crystal plasticity (CP) for the ferrite phase with J2 isotropic plasticity for the martensite phase. The effects of transformation-induced geometrically necessary dislocations (GNDs) are incorporated through an additional thermal strain field representing the volumetric expansion associated with martensitic transformation. Phase-field damage evolution, coupled with the modified Mohr–Coulomb (MMC) failure criterion, is implemented through the HETVAL user subroutine to account for stress triaxiality- and Lode-parameter-dependent fracture behavior.
mod_alloys.f defines and updates the crystallographic and material parameters of the ferrite phase used in the crystal plasticity model. 

umat_matj2.f - a user-defined material subroutine (UMAT) implementing J2 isotropic plasticity for both ferrite and martensite phases in dual-phase (DP) steels to facilitate the calibration of crystal plasticity (CP) parameters for ferrite against experimental mechanical responses.

# Data: 
TensileRaw.xlsx - the experimental tensile test data for dual-phase (DP) steels with fibrous and coarse martensite morphologies. These datasets are used to calibrate and validate the corresponding finite element simulations by comparing experimentally measured and numerically predicted stress–strain responses, ensuring that the constitutive models accurately reproduce the macroscopic mechanical behavior of both microstructures.

# Inp:
The Abaqus input files (.inp) generated from electron backscatter diffraction (EBSD) data for dual-phase (DP) steels with fibrous and coarse martensite morphologies. These files are used to perform finite element simulations under different ferritic crystallographic texture configurations.

