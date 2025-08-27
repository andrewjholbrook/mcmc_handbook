using Pkg;
pkg"activate .";
#PKG_ROOT=ENV["PKG_ROOT"];
#pkg"activate $PKG_ROOT";

#scen="154835"; computeGradients=true;
#scen="scal_var"; computeGradients=false;
scen="consistency"; computeGradients=true;

println("Starting setup...");
include("setup_$(scen).jl");
println("Finished setup.");

s0 = sTrue[sampInd];
#s0 = zeros(nSampInd);

#tic();
#samples,obs,lpdfs,ar = adRunMcmc(mcmcP, s0; computeGradients=computeGradients, verbose=2);
#mcmcTime = toq();
mcmcTime = @elapsed samples,obs,lpdfs,ar = adRunMcmc(mcmcP, s0; computeGradients=computeGradients, verbose=2);

samples=permutedims(hcat(samples...),[2,1]);
obs=permutedims(hcat(obs...),[2,1]);

#h5fl="$(scen).h5";
cnt=1;
h5fl="$(scen)_$(cnt).h5";
while (isfile(h5fl))
  global cnt, h5fl;
  h5fl="$(scen)_$(cnt).h5";
  cnt+=1;
end

println("Writing output to $(h5fl)...");
h5write(h5fl,"samples",samples);
h5write(h5fl,"obs",obs);
h5write(h5fl,"lpdfs",lpdfs);
h5write(h5fl,"ar",ar);

acceptRatio = sum(ar[1:mcmcP.nsamp])/mcmcP.nsamp;
h5write(h5fl,"acceptRatio",acceptRatio);

sampPerSec = mcmcP.nsamp / mcmcTime;
h5write(h5fl,"mcmcTime",mcmcTime);
h5write(h5fl,"sampPerSec",sampPerSec);
println("Wrote $(h5fl).");
