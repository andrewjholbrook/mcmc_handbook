using Printf, HDF5, StatsBase
using DataFrames, CSV, StatsPlots
using DataFramesMeta
using Measures

#fix errors for headless plotting
#GKS: can't connect to GKS socket application
ENV["GKSwstype"] = "100";
include("../../plot/plotSave.jl");

rawCsvFile = "AC_AR_raw.csv";
csvFile    = "AC_AR.csv";
if isfile(rawCsvFile)
  println("Reading data from raw csv file: $(rawCsvFile)");
  df = DataFrame(CSV.File(rawCsvFile));
else
  files = [ 
              [ @sprintf("/projects/SIAllocation/lolmcGalerkin/lolmcGalerkin_%03d.h5",i) for i=[13:14;17:42] ];
              [ @sprintf("/projects/SIAllocation/2020juqex2/2020juqex2_%03d.h5",i) for i=12:49 ];
              [ @sprintf("/projects/SIAllocation/lolmcNN/lolmcNN_%03d.h5",i) for i=5:48 ]
          ];
  
  mcmcs  = [ "" for i=1:length(files) ];
  gnorms = zeros(Int,length(files));
  ars    = zeros(length(files));
  acs    = zeros(length(files));
  msjds  = zeros(length(files));
  success = [ false for i=1:length(files) ];
  
  println("Reading data from files...");
  for i = 1:length(files)
      file = files[i];
      if isfile(file)
          try
              gnorms[i] = h5read(file,"maxnormGalerkin");
          catch e; end
  
          try
              mcmcs[i]  = h5read(file,"mcmc"); 
              ars[i]    = h5read(file,"acceptRatio"); 
              lpdfs     = h5read(file,"lpdfs"); 
              samples   = h5read(file,"samples");
              nsamp     = h5read(file,"sampComplete");
              
              llhAC = autocor(lpdfs[:,2]);
              meanJumpDist = mean(sum((samples[2:end,:]-samples[1:end-1,:]).^2,dims=2));
  
              acs[i]    = llhAC[30];
              msjds[i]  = meanJumpDist;
  
              @printf("    %20s %16s %2d %8d %8.4f %8.4f %8.4f\n",basename(file),mcmcs[i],gnorms[i],nsamp,ars[i],acs[i],msjds[i]);
              success[i] = true;
          catch e
              println("    Error getting data from $(file).");
              println("    File may be malformed or incomplete.");
              println("    Error was: $(e)");
          end
  
      else
          println("    File not found: $(file). Skipping...");
      end
  end

  println("Creating initial dataframe...");
  df = DataFrame((File=[ basename(file) for file in files ],MCMC=mcmcs,Gnorm=gnorms,AR=ars,AC30=acs,MSJD=msjds));

  #remove files that failed to load
  df = df[success,:];
  
  println("Writing initial csv...");
  CSV.write(rawCsvFile,df);
  println("Wrote: $(rawCsvFile)");
end

println("Processing...");
#df.Scenario  = [ Meta.parse(split(str,"_")[1]) for str in df.File ];
df.Scenario  = [ split(str,"_")[1] for str in df.File ];

#average AR results by mcmc column
df = groupby(df,[:MCMC,:Gnorm,:Scenario]);
df = combine(df,:AR => mean, :AC30 => mean, :MSJD => mean);

#split mcmc column
df.Method  = [ Meta.parse(split(str,"|")[1]) for str in df.MCMC ];

df.Eps .= 0.0; #missing;
df.I   .= 0; #missing;
df.Tau .= 0.0; #missing;
df.H   .= 0.0; #missing;

hmcRows  = ((df.Method .== :hmc) .| (df.Method .== :lolhmc));
malaRows = (df.Method .== :mala);
df.Eps[hmcRows] = [ eval(Meta.parse(split(str,"|")[2])) for str in df.MCMC[hmcRows] ];
df.I[hmcRows]   = [ eval(Meta.parse(split(str,"|")[3])) for str in df.MCMC[hmcRows] ];
df.Tau[hmcRows] = df.Eps[hmcRows] .* df.I[hmcRows];
df.H[malaRows]  = [ eval(Meta.parse(split(str,"|")[2])) for str in df.MCMC[malaRows] ];

df.AC30_norm .= 0.0; #missing;
df.MSJD_norm .= 0.0; #missing;

df.AC30_norm[hmcRows] = df.AC30_mean[hmcRows] ./ df.I[hmcRows];
df.MSJD_norm[hmcRows] = df.MSJD_mean[hmcRows] ./ df.I[hmcRows];

display(df)

println("Writing processed csv...");
CSV.write(csvFile,df);
println("Wrote: $(csvFile)");

#function triplePlot(pltDf::DataFrame, xCol::Symbol, xLab::String, grpCol::Symbol, grpLab::String; size=(1200,400))
function triplePlot(pltDf::DataFrame, xCol::Symbol, xLab::String, grpCol, grpLab; size=(1200,400), yCols=[:AR_mean,:AC30_mean,:MSJD_mean], yLabs=["Acceptance Ratio","Autocorrelation in \\Phi","Mean Squared Jumping Distance"])
    p = plot(layout=(1,3),size=size,margin=5mm)#,leg=:outerbottom)
    yCol = yCols[1];
    @with pltDf scatter!(p[1],
        $xCol,
        $yCol,
        group = $grpCol,
        ylim=(0,1),
        xlab=xLab,
        ylab=yLabs[1]
    )
    yCol = yCols[2];
    @with pltDf scatter!(p[2],
        $xCol,
        $yCol,
        group = $grpCol,
        ylim=(0,1),
        xlab=xLab,
        ylab=yLabs[2]
    )
    yCol = yCols[3];
    @with pltDf scatter!(p[3],
        $xCol,
        $yCol,
        group = $grpCol,
        xlab=xLab,
        ylab=yLabs[3]
    )
    for i=1:length(p)
      if grpLab != ""
        plot!(p[i],legendtitle=grpLab);
      else
        plot!(p[i],legend=false);
      end 
    end
    return(p)
end

p = triplePlot(df[(df.Method.==:hmc   ).&(df.Scenario.=="2020juqex2"   ),:], :Tau, "End Time (\\tau)", :I, "Integration Steps");
plotSave(p,"figures/hmc_3plots");

p = triplePlot(df[(df.Method.==:hmc   ).&(df.Scenario.=="2020juqex2"   ),:], :Tau, "End Time (\\tau)", :I, "Integration Steps"; yCols=[:AR_mean,:AC30_norm,:MSJD_norm]);
plotSave(p,"figures/hmc_3plots_norm");

p = triplePlot(df[(df.Method.==:hmc   ).&(df.Scenario.=="lolmcNN"      ),:], :Tau, "End Time (\\tau)", :I, "Integration Steps");
plotSave(p,"figures/nn_lolhmc_3plots");

p = triplePlot(df[(df.Method.==:lolhmc).&(df.Scenario.=="lolmcGalerkin"),:], :Tau, "End Time (\\tau)", :Gnorm, "Galerkin Norm");
plotSave(p,"figures/gal_lolhmc_3plots");

p = triplePlot(df[(df.Method.==:mala  ).&(df.Scenario.=="lolmcGalerkin"),:], :H, "Step Size (h)", :Gnorm, "Galerkin Norm");
plotSave(p,"figures/gal_lolmala_3plots");

p = triplePlot(df[(df.Method.==:mala  ).&(df.Scenario.=="lolmcNN"      ),:], :H, "Step Size (h)", :Scenario, ""); #hack to make no labels
plotSave(p,"figures/nn_lolmala_3plots");

