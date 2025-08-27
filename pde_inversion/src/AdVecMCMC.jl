module AdVecMCMC

using LinearAlgebra
using Printf
using SparseArrays
using SpectralDiscrete2D
using InfDimMCMC
using AdvectionDiffusion

export adPointData;
export adSpectralData;
export adScalVarData;

export adSpectralGradPhi;
export adFemGradPhi;

export kraichnanEnergy;

export adRunMcmc;

include("adPointData.jl");
include("adSpectralData.jl");
include("adScalVarData.jl");
include("adSpectralGradPhi.jl");
include("adFemGradPhi.jl");
include("kraichnanEnergy.jl");
include("adRunMcmc.jl");

end # module
