function adRunMcmc(mcmcP::mcmcProb, cur::mcmcSample; computeGradients=true, verbose=0)

  #initialize
  mcmcFillSample(cur,mcmcP;computeGradients=computeGradients);
  
  #run mcmc forward
  samples = Array{typeof(cur.samp)}(undef,mcmcP.nsamp);
  obs     = Array{typeof(cur.obs)}(undef,mcmcP.nsamp);
  ar      = Array{UInt8}(undef,mcmcP.nsamp);
  lpdfs   = Array{Float64}(undef,mcmcP.nsamp,3);

  #burnin
  (verbose>0) && @printf("\n-----------------------Starting Burn In-----------------------\n");
  for i=1:mcmcP.nburn
    (verbose>0) && @printf("\nBurn in Step #%d:\n",i);
    cur,_,_ = mcmcP.step(cur, mcmcP; verbose=verbose);
  end

  #(verbose>0) && @printf "\n------------Starting Metropolis-Hastings Sampling------------\n"
  (verbose>0) && @printf("\n----------------------Starting Sampling-----------------------\n");
  for i=1:mcmcP.nsamp
    (verbose>0) && @printf("\nSample #%d:\n",i);
    cur,ar[i],can = mcmcP.step(cur, mcmcP; verbose=verbose);
    samples[i] = cur.samp;
    obs[i]     = cur.obs;
    lpdfs[i,1] = cur.prLogPdf;
    lpdfs[i,2] = cur.llLogPdf;
    lpdfs[i,3] = cur.postLogPdf;
  end
  (verbose>0) && @printf("\n------------Sampling Complete------------\n");

  return samples, obs, lpdfs, ar;
end

function adRunMcmc(mcmcP::mcmcProb, s0=rand(mcmcP.prior); computeGradients=true, verbose=0)
  #initialize
  cur = mcmcSample();
  cur.samp = s0;

  return adRunMcmc(mcmcP, cur; computeGradients=computeGradients, verbose=verbose)
end
