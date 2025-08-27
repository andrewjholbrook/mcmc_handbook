#generate plots for the scenario
using HDF5
using InfDimMCMC
using Printf
using LinearAlgebra
#using StatsBase
#using Statistics
#using Distributions
#using AdVecMCMC

#if the scenario hasn't been run (use lpdfs as a proxy for this), load it
#also allow re-downloading of the data by specifying rerun=true
sthin=1;
if ((!@isdefined lpdfs) || ((@isdefined rerun) && rerun))
  samples, _, lpdfs, ar, files = assembleChain(outFile;verbose=true,sampleThin=sthin,sampleCols=1:8,obsCols=1:2);
  #ENV["ADR_ROOT"]="/projects/SIAllocation";
  #include("setup.jl");
  # f = h5open(outFile);
  # #maxnorm  = read(f,"maxnorm");
  # #maxnormv = read(f,"maxnormV");
  # #nsstd    = read(f,"noiseStd");
  # lpdfs    = read(f,"lpdfs");
  # samples  = read(f,"samples");
  # close(f);
end

using Plots
using Plots.Measures

include("../../plot/adComputeTotVarEvolve.jl");
include("../../plot/adPlotTotVarEvolve.jl");
include("../../plot/adPlotHistMatrix.jl");
include("../../plot/plotSave.jl");

#fix errors for headless plotting
#GKS: can't connect to GKS socket application
ENV["GKSwstype"] = "100"

#read true samples (make sure to truncate the empty first line)
samplesTrue = h5read("/projects/SIAllocation/runs/206380/adr_bayes.h5","samples")[:,2:end];
#convert to new format
samplesTrue[:,1:2:end] .*=  sqrt(2);
samplesTrue[:,2:2:end] .*= -sqrt(2);

gr();
#pltScl = 2.0;

#autocorrelation
llhAC = autocor(lpdfs[:,2]);
vACidx = 1:4; vAC = [ autocor(samples[:,i]) for i=vACidx ];
p = plot(ylim=(0,1), xlab = "Lag", ylab = "Autocorrelation",size=(450,300));
plot!(p,llhAC,lab="llh");
for i=1:length(vAC); plot!(p,vAC[i],lab="v$(vACidx[i])"); end
plotSave(p,replace(outFile,".h5"=>"_ac"));

#trace plot
nthin=1000; idx=1:nthin:size(samples,1);
gr(); 
p = plot(layout=(4),xlab="sample",legend=false); 
for i=1:4
    plot!(p[i],idx,samples[idx,i],ylab="v$(i)"); 
end
plotSave(p,replace(outFile,".h5"=>"_trace"));

#histogram vs. truth
nthin=1;
bins=-3:0.1:3;
midpts = 0.5*( bins[2:end] + bins[1:end-1] );
gr(); 
p = plot(layout=(4),xlab="sample",legend=false); 
for i=1:4
#     stephist!(p[i],samples[1:nthin:end,i],ylab="v$(i)",bins=100,normalize=:pdf); 
#     stephist!(p[i],samplesTrue[:,i],ylab="v$(i)",bins=100,normalize=:pdf);
    h1 = fit(Histogram, samples[1:nthin:end,i], bins);
    h1 = normalize(h1, mode=:pdf);
    h2 = fit(Histogram, samplesTrue[:,i], bins);
    h2 = normalize(h2, mode=:pdf);
    plot!(p[i],midpts,h1.weights);
    plot!(p[i],midpts,h2.weights);
#     @printf("TV %d: %10.6f\\n",i,0.5*dot(bins[2:end] - bins[1:end-1],abs.(h2.weights-h1.weights)))
end
plotSave(p,replace(outFile,".h5"=>"_hist"));

#convergence in total variation
nsteps = 100;
tvFile = replace(outFile,".h5"=>"_totVarEvolve.h5");
if !isfile(tvFile)
    adComputeTotVarEvolve( samples, samplesTrue; outFile=tvFile, comp=1:4, nsteps=nsteps);
end
adPlotTotVarEvolve(tvFile,replace(outFile,".h5"=>"_totvar"));

#2d histogram
nthin=1000;
# adPlotHistMatrix(samples[1:nthin:end,1:4],"hist2d_lolmc";nbins=100,rng=[-3 3],comp=1:4,size=(800,800),
#     margin=0.0mm,xticks=-3:3,yticks=-3:3, label=["v$(i)" for i=1:4]);
bins=[ range(-3,stop=3,length=100) range(-3,stop=3,length=100) ];
p = histmatrix(samples[1:nthin:end,1:4];bins=bins,size=(600,600),
    margin=0.0mm,xticks=-3:3,yticks=-3:3, label=["v$(i)" for i=1:4]);
plotSave(p,replace(outFile,".h5"=>"_histmatrix"));


