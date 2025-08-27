using Pkg;
pkg"activate .";

using HDF5


ADR_ROOT=ENV["ADR_ROOT"];

oldDataFile = "/home/jkrometi/math/adr_bayes/julia/data/point_twohump_012.h5";
newDataFile = ADR_ROOT*"/data/point_twohump_012.h5";

#read these values from the old datafile
maxnorm = h5read(oldDataFile,"knorm"); 
kappa   = h5read(oldDataFile,"kappa"); 
t0      = h5read(oldDataFile,"t0");
dt      = h5read(oldDataFile,"dt");
tf      = h5read(oldDataFile,"tf"); 
dataTX  = h5read(oldDataFile,"dataTX");
y       = h5read(oldDataFile,"y");

nk      = length(h5read(oldDataFile,"vk_true"))+1;


#need to set these up manually due to a change in how we're representing them
th0     = zeros(nk); th0[1] = 0.5; th0[3:2:5] .= -1/(4*sqrt(2));
sTrue   = zeros(198); sTrue[3] = -4*sqrt(2); sTrue[5] = 4*sqrt(2); sTrue = sTrue[3:end];


#write everything
f = h5open(newDataFile,"w");
write(f,"maxnorm",maxnorm);
write(f,"kappa",kappa);
write(f,"t0",t0);
write(f,"dt",dt);
write(f,"tf",tf);
write(f,"dataTX",dataTX);
write(f,"y",y);
write(f,"th0",th0);
write(f,"sTrue",sTrue);
close(f);

println("Wrote: $newDataFile");

