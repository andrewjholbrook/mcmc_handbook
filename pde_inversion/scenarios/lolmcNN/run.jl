#using Pkg;
#pkg"activate .";
#PKG_ROOT=ENV["PKG_ROOT"];
#pkg"activate $PKG_ROOT";
using HDF5

#Parse command line arguments
using ArgParse
apSettings = ArgParseSettings();
@add_arg_table apSettings  begin
  "--scen"
    help = "scenario to run"
    required = false
  "--restartfile"
    help = "file to restart from"
    required = false
  "--init"
    help = "initialization for chain (ignored if --restartfile is specified)"
    required = false
    default = "mle"
  "--mcmc"
    help = "mcmc method and parameters (comma-delimited)"
    required = false
  "--nsamp"
    help = "number of samples to run"
    arg_type = Int
    required = false
  "--nburn"
    help = "number of burnin samples to run"
    arg_type = Int
    required = false
  "--modelfile"
    help = "model file"
    required = false
end
args = parse_args(apSettings,as_symbols=true);

#Print command line arguments and convert them to global variables
#(This allows the setup scripts to be run without having to pass in a dictionary)
println("Parsed args:")
for (arg,val) in args
  if val != nothing
    println("  $arg  =>  $val"); #print
    @eval (($arg) = ($val));     #global variable
  end
end
if (@isdefined restartfile)
  modelfileTmp = h5read(restartfile,"modelfile");
  if (@isdefined modelfile) && (modelfileTmp != modelfile)
    error("modelfile is specified ($(modelfile)) but is being overwritten by modelfile from restartfile ($(modelfileTmp))!");
  end
  modelfile = modelfileTmp;
  if (!@isdefined mcmc)
    mcmc = h5read(restartfile,"mcmc");
  end
end


println("Starting setup...");
include("setup.jl");
println("Finished setup.");


outDir="$(ADR_ROOT)/$(scen)";
if ( ! isdir(outDir) ) 
  println("Output directory $(outDir) does not exist. Creating...");
  mkdir(outDir);
end


## Setup output file ##
mkfilename(cnt;outDir=outDir,scen=scen) = @sprintf("%s/%s_%03d.h5",outDir,scen,cnt);
cnt=1;
outFile=mkfilename(cnt);
while (isfile(outFile))
  global cnt, outFile;
  outFile=mkfilename(cnt);
  cnt+=1;
end
println("Writing output to $(outFile)...");
h5write(outFile,"datafile",datafile);
h5write(outFile,"modelfile",modelfile);
h5write(outFile,"dataDim",data.dim);
(@isdefined nburn) && h5write(outFile,"nburn",nburn);
(@isdefined nsamp) && h5write(outFile,"nsamp",nsamp);
(@isdefined mcmc ) && h5write(outFile,"mcmc",mcmc);
#save command line arguments
for (key,val) in args
  if val != nothing
    h5write(outFile,"args/$(key)",val);
  end
end


#Initial sample (use method specified by --init unless --restartfile is provided)
function restartSample(filename)
  #return h5read(filename,"samples")[end,:];
  f = h5open(filename);
  dset = f["samples"];
  s = dset[end,:][:];
  close(f);
  return s;
end
if (@isdefined restartfile)
  s0 = restartSample(restartfile);
  h5write(outFile,"restartfile",restartfile);
  h5write(outFile,"init","restart");
  println("initializing chain from $(restartfile)");
else
  if (init == "mle")
    trainFile="/projects/SIAllocation/lolmc/juqex2_traindata.h5";
    trSamples = h5read(trainFile,"samples");
    trPots = h5read(trainFile,"pots");
    s0 = trSamples[argmin(trPots),3:end]; 
  elseif (init == "rand")
    s0 = rand(mcmcP.prior);
  elseif (init == "zero")
    s0 = zeros(nSampInd);
  elseif (init == "true")
    s0 = sTrue[sampInd]; 
  end
  h5write(outFile,"init",init);
  println("initializing chain using method $(init)");
end



## Run ##
#mcmcTime = @elapsed samples,obs,lpdfs,ar = adRunMcmc(mcmcP, s0; computeGradients=computeGradients, verbose=2);
mcmcTime = @elapsed mcmcRun(mcmcP, s0; verbose=3, outFile=outFile);



## Post process ##

ar      = h5read(outFile,"ar");
lpdfs   = h5read(outFile,"lpdfs");

acceptRatio = sum(ar[1:mcmcP.nsamp])/mcmcP.nsamp;
sampPerSec = mcmcP.nsamp / mcmcTime;
secPerSamp = mcmcTime / mcmcP.nsamp;

h5write(outFile,"acceptRatio",acceptRatio);
h5write(outFile,"mcmcTime",mcmcTime);
h5write(outFile,"sampPerSec",sampPerSec);
h5write(outFile,"secPerSamp",secPerSamp);

println("Wrote $(outFile).");

@printf("MCMC Time for %d samples: %.2f seconds (%.4f seconds/sample)\n",mcmcP.nsamp,mcmcTime,secPerSamp);
@printf("Acceptance ratio: %.4f\n\n",acceptRatio);

stp=max(1,round(Int,mcmcP.nsamp/20));

@printf("\nSummary of run:\n");
@printf("%9s  %9s: %8s ", "start", "end", "ar (avg)"); 
@printf("%16s %16s %16s\n", "logprior", "loglikelihood", "logposterior");
for i=0:stp:mcmcP.nsamp-1
  @printf("%9d -%9d:",i+1,i+stp);
  @printf(" %8.4f",mean( ar[i+1:i+stp] ));
  @printf(" %16.6f %16.6f %16.6f\n", lpdfs[i+stp,1], lpdfs[i+stp,2], lpdfs[i+stp,3]);
end

# ## Plots ##
# include("plot.jl");
