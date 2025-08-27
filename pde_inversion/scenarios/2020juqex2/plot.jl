using Plots
using Plots.Measures

include("../../plot/plotSave.jl");

#fix errors for headless plotting
#GKS: can't connect to GKS socket application
ENV["GKSwstype"] = "100"

gr();

hmgusz=32;
hmtksz=24;

include("../../plot/adPlotHistMatrix.jl");
include("../../plot/adPlotVTrace.jl");

adPlotHistMatrix(outFile;outFile=replace(outFile,".h5"=>"_histmatrix"),nbins=100,rng=[-2 2],comp=1:8,size=(1600,1600),guidefontsize=hmgusz,tickfontsize=hmtksz,margin=0.0mm,xticks=-1:1:1,yticks=-1:1:1);

adPlotVTrace(outFile;outFile=replace(outFile,".h5"=>"_vtrace"),comp=1:4);
