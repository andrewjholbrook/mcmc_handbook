using Pkg;
pkg"activate .";

include("setup.jl");

mcmcP.computeGradients = true;
s = mcmcSample();
s.samp = rand(mcmcP.prior);
mcmcFillSample(s,mcmcP);


sTmp = mcmcSample();
sTmp.samp = copy(s.samp);
mcmcP.computeGradients = false;

del = 0.0001;
ncomp = 10;
comp = sort(rand(1:length(s.samp),ncomp));
@printf("%6s: %16s %16s %16s\n","comp","adj","findiff","adj/fd");
for ci=comp
  sTmp.samp[ci] = s.samp[ci] + del;
  mcmcFillSample(sTmp,mcmcP);
  potp = sTmp.pot;
  
  sTmp.samp[ci] = s.samp[ci] - del;
  mcmcFillSample(sTmp,mcmcP);
  potn = sTmp.pot;
  #@printf("%6d: %15.8f %15.8f %15.8f\n",ci,s.pot,potp,potn);

  gp = (potp - potn)/(2*del);
  @printf("%6d: %16.8f %16.8f %16.8f\n",ci,s.gradPot[ci],gp,s.gradPot[ci]/gp);
end

