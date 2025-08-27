using LinearAlgebra
using Printf
using Distributions
using HDF5

using CUDA
using Flux
#using Flux: @epochs, train!
using BSON: @save, @load
#using ProgressMeter
using StatsBase
using Statistics

using SpectralDiscrete2D
using AdvectionDiffusion
using InfDimMCMC
using AdVecMCMC

ADR_ROOT=ENV["ADR_ROOT"];

#defaults (typically overwritten by arguments to run.jl)
def_datafile  =ADR_ROOT*"/data/point_twohump_012.h5";
def_mcmc  = "lolhmc|2^-3|2";
def_nsamp = 10;
def_nburn = 0;
def_modelfile = ADR_ROOT*"/lolmc/trained_models/model_001.bson";

## ADVECTION-DIFFUSION PROBLEM DEFINITION: KAPPA, INITIAL CONDITION, OBSERVATIONS/DATA ##
ad = adProb();

datafile = (@isdefined datafile) ? datafile : def_datafile;
@printf("Reading problem setup from datafile: %s\n",datafile);
maxnorm   = 8;#h5read(datafile,"maxnorm"); 
ad.kappa  = h5read(datafile,"kappa"); 
ad.t0     = h5read(datafile,"t0");
ad.dt     = h5read(datafile,"dt");
ad.tf     = h5read(datafile,"tf"); 
ad.th0    = h5read(datafile,"th0"); 
dataTX    = h5read(datafile,"dataTX");
y         = h5read(datafile,"y");
#nsStd     = h5read(datafile,"nsStd");

#ad.dt = 0.001;
ad.t      = collect( ad.t0:ad.dt:ad.tf ); ad.nt = length(ad.t);

#kdisc = specDiscrete2D(maxnorm;multiplicity=true,computeKdiff=true);
kdisc = specDiscrete2D(maxnorm;multiplicity=true,kDiffFile=@sprintf("%s/lookup/kdiff_kd%03d.h5",ADR_ROOT,maxnorm));

datafileMaxNorm = h5read(datafile,"maxnorm");
if datafileMaxNorm > kdisc.maxnorm
  @printf("datafile maxnorm > kdisc.maxnorm (%d < %d). Truncating th0 to match...\n", datafileMaxNorm, kdisc.maxnorm);
  ad.th0    = ad.th0[1:kdisc.nk];
elseif datafileMaxNorm < kdisc.maxnorm
  @printf("datafile maxnorm < kdisc.maxnorm (%d < %d). Padding th0 with zeros to match...\n", datafileMaxNorm, kdisc.maxnorm);
  ad.th0    = [ad.th0; zeros(kdisc.nk - length(ad.th0))];
end


data=adPointData(ad,kdisc,dataTX,y;computeL2Kern=false);

## Setup NN ##
modelfile = (@isdefined modelfile) ? modelfile : def_modelfile;
@load modelfile model;
model = model |> gpu;


## Setup MCMC Problem ##

mcmcP = mcmcProb();

# Number of samples # 
mcmcP.nburn = (@isdefined nburn) ? nburn : def_nburn;
mcmcP.nsamp = (@isdefined nsamp) ? nsamp : def_nsamp;

# Define sampler # 
mcmc = (@isdefined mcmc) ? mcmc : def_mcmc;
mcmcSetSampler(mcmcP,mcmc);

## Setup sample space ##
maxnormV = maxnorm;
kdiscV = specDiscrete2D(maxnormV;multiplicity=true,kDiffFile=@sprintf("%s/lookup/kdiff_kd%03d.h5",ADR_ROOT,maxnormV));
sampInd   = 3:kdiscV.nk;   #Components to consider
sampNoInd = setdiff(1:kdiscV.nk,sampInd);   #Known k's
nSampInd  = length(sampInd);

# Prior #
prStd = sqrt.(kraichnanEnergy(kdisc.kn));
prStd[kdisc.kn .== 0.0] .= 0.25;
prStd = prStd[sampInd];
prStd .*= sqrt(2); #conversion to match prior from 2020juqex2
mcmcP.prior = MvNormal( zeros(nSampInd), prStd );#sqrt.(kraichnanEnergy(kdisc.kn)[sampInd]) );

# Map from samples to vector components #
let nk=kdiscV.nk, sampInd=sampInd
  function padZeros(s)
    p = zeros(nk);
    p[sampInd] = s;
    return p;
  end
  InfDimMCMC.mcmcSampToParamMap(s) = padZeros(s.samp);
end
let nk=kdiscV.nk
  gspMat = [ zeros(2,nSampInd); I ];
  InfDimMCMC.mcmcGradSampToParamMap(s) = gspMat;
end

#make a v with all non-zero components and compute the template
ad.v = zeros(maximum(sampInd)); 
Ip,Tp,Inn,Tn,T0 = adSpectralMatrixTemplate(ad,kdisc);

# Forward map and observations #
let Ip=Ip, Tp=Tp, Inn=Inn, Tn=Tn, T0=T0
  function adSolve(v)
    ad.v = v;
    return adSpectral(ad,Ip,Tp,Inn,Tn,T0);
  end
  InfDimMCMC.mcmcForwardMap(s) = adSolve(s.param);
end

# Observation map #
let kdisc=kdisc, data=data
  local specBaseXY = spectralBasis(kdisc,data.dataTX[:,2:3]);

  function obs(thk)
    return diag( data.tmat * (thk*specBaseXY') );
  end
  InfDimMCMC.mcmcObsMap(s) = obs(s.sol);
end

# Potential map #
let ad=ad, kdisc=kdisc, data=data, Ip=Ip, Tp=Tp, Inn=Inn, Tn=Tn, T0=T0
  nsStd = 0.125; #honestly unclear if this is std or variance
  noiseMean = zeros(data.dim);
  Q0 = MvNormal(noiseMean,nsStd);
  # #dPhidG(Gj) = -invcov(Q0)*(data.y-Gj);
  # dPhidG(Gj) = -Q0.Σ.inv_value*(data.y-Gj); #faster for simple covariance structure
  # 
  # #we'll need this to compute the adjoint forcing
  # dk = deltak(kdisc,data.dataTX[:,2:3]) ./ ad.dt;

  #mcmcP.potMap = ( s -> ( -logpdf(Q0,data.y-s.obs) + logpdf(Q0,data.y) ) );
  potMap = ( s -> ( -logpdf(Q0,data.y-s.obs) ) );
  InfDimMCMC.mcmcPotMap(s) = potMap(s);
end

# Gradient of potential map #
let model=model

  function gradPot(v)
    #return model(v |> gpu) |> cpu;
    return model(reshape(v,length(v),1) |> gpu)[:] |> cpu;
  end

  InfDimMCMC.mcmcGradPotMap(s) = gradPot(s.samp);
end

# #TEST
# s = mcmcSample(); s.samp = rand(mcmcP.prior); mcmcFillSample(s,mcmcP);
# #@printf("norm(G(sTrue)-y) = %12.8f\n",norm(s.obs-y));
