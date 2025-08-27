using HDF5, Distributions;

include("specDiscrete2D.jl");    #2D spectral discretization
include("adProb.jl");            #adProb definition
include("vectorField.jl");       #vector field functions
include("ode.jl");               #ode simulators
include("adSpectral.jl");        #spectral AD solver
include("timeInterp.jl");        #time interpolation
include("util.jl");              #utilities
include("kraichnanEnergy.jl");   #kraichnan energy spectrum

## SETUP ##

ADR_ROOT=ENV["ADR_ROOT"];

if (!isdefined(:h5fl))
  h5fl="$ADR_ROOT/data/sparse.h5";
  @printf("h5fl not defined. writing to default of %s...\n", h5fl);
end

if (!isdefined(:method))
  method = "numerical";
  @printf("method not defined. using default value of %s...\n", method);
end

if method=="numerical"

  #SPECTRAL DISCRETIZATION
  if (!isdefined(:kdisc))
    @printf("kdisc not defined. using default value for knorm...\n");
    knorm = 32;           #affects how many k's to use
    kdisc = specDiscrete2D(knorm);
    
    #get map k_i-k_j -> k (read from file)
    kdfl=@sprintf("%s/lookup/kdiff_kd%03d.h5",ADR_ROOT,knorm);
    if (!isfile(kdfl))
      error("ERROR: kdiff file $kdfl does not exist\n");
    end
    kdisc.kdiff = h5read("$kdfl","kd");
  end
  @printf("using knorm=%d (nk=%d)...\n",kdisc.knorm,kdisc.nk);

  if (!isdefined(:dataTX))
    dataDim=4096;
    #dataTX = [0.5 0.5 0.5; 0.5 0.5 1.0; 1.0 0.5 0.5; 1.0 0.5 1.0];
    dataTX = rand(dataDim,3);
    dataTX[:,1] = sort(dataTX[:,1]); #order by time for convenience
    @printf("dataTX not defined. randomly generating:\n");
    #show(dataTX);
    for i=1:size(dataTX,1); @printf("%8.6f  %8.6f  %8.6f\n", dataTX[i,1], dataTX[i,2], dataTX[i,3]); end
  end
  if (!isdefined(:kappa))
    #kappa = 3e-5; #from Run 1, Table I, chen1998simulations
    kappa = 0.282; #from Wikipedia for Water in Air
    @printf("kappa not defined. using default value...\n");
  end

  #THETA0
  if (!isdefined(:th0k))
    #th0k = zeros(Complex{Float64},kdisc.nk); th0k[1]=0.5; th0k[2:3]=0.25;    #0.5*cos(2*pi*y)+0.5
    th0k = zeros(Complex{Float64},kdisc.nk); th0k[1]=0.5; th0k[2:5]=-0.125;    #0.5-0.25*cos(2*pi*x)-0.25*cos(2*pi*y)
    @printf("th0k not defined. using default value...\n");
  end

  #VK_TRUE
  if (!isdefined(:vk_true))
    #@printf("vk_true not defined. using default value...\n");
    #vk_true = zeros(Complex{Float64},kdisc.nk); vk_true[2:2:4]=1; vk_true[3:2:5]=-1;
    #vk_true = zeros(Complex{Float64},kdisc.nk); vk_true[2]=-1; vk_true[4]=1;  #from original quad hump distribution
    #vk_true = zeros(Complex{Float64},kdisc.nk); vk_true[2]=-4; vk_true[4]=4;
    @printf("vk_true not defined. drawing at random from the Kraichnan prior...\n");
    mu0_ind  = 2:kdisc.nk;
    mu0_mean = zeros(kdisc.nk);
    mu0_var  = diagm( kraichnanEnergy(kdisc.kn[mu0_ind]) );
    mu0      = MvNormal(mu0_mean[mu0_ind],mu0_var); #mu0 is multivariate normal
    vk_true = zeros(Float64,kdisc.nk);
    vk_true[mu0_ind] = rand(mu0);
    vk_true = rl2vk(vk_true);
  end

  #Build AD problem
  ad = adProb();
  ad.kappa = kappa;
  ad.t0 = 0.0; ad.dt=0.01; ad.tf=1.0;
  ad.t=collect( 0:ad.dt:ad.tf ); ad.nt = length(ad.t);
  ad.th0 = th0k;
  ad.v = vk_true;

  #print the adProb
  adPrint(ad);

  thk = adSpectral(ad,kdisc);
  
  #spatial bases evaluated at dataTX points
  specBaseXY = e2piikx(kdisc.k,dataTX[:,2:3]);
  #time interpolation matrix
  tmat = timeInterp(ad.t,dataTX[:,1]);

  G_true = zeros( dataDim );
  for i=1:size(tmat,1); 
    G_true[i] = ( tmat[i,:].'*real(thk*specBaseXY[i,:]) )[1]; 
  end

end

if method=="exact1"
  #exact solution for:
  #  th0   = cos(2*pi*y)
  #  v     = [0, v_y]
  th0k = zeros(kdisc.nk); th0k[1]=0.5; th0k[2:3]=0.25;    #0.5*cos(2*pi*y)+0.5

  if (!isdefined(:v))
    v = [0 1];
    @printf("v not defined. using default value: [%8.6f, %8.6f]\n", v[1], v[2]);
  end
  if (!isdefined(:kappa))
    kappa = 0.025;
    @printf("kappa not defined. using default value of %8.6f...\n", kappa);
  end
  
  k = 1;
  thexf(t,x) = 0.5+0.5*exp(-4*pi^2 * kappa^2 * k^2 * t ) .* cos( 2*pi*k*(x[:,2] - v[2]*t) );

  G_true = thexf(dataTX[:,1],dataTX[:,2:3]);

  DESCRIPTION="Exact solution for th0=cos(2*pi*y); v=[0,v_y]";
  HOWTO="v=[0 1]; k=1; thexf(t,x) = 0.5+0.5*exp(-4*pi^2 * kappa^2 * k^2 * t ) .* cos( 2*pi*k*(x[:,2] - v[2]*t) )";
end

if method=="exact2"
  #  th0   = cos(2*pi*y)
  #  kappa = 0
  #  v     = [0 sin(2*pi*x)];
  kappa = 0;
  th0k = zeros(kdisc.nk); th0k[1]=0.5; th0k[2:3]=0.25;    #0.5*cos(2*pi*y)+0.5
  thexf(t,x) = 0.5+0.5*cos( 2*pi*(x[:,2]-sin(2*pi*x[:,1]).*t) );

  G_true = thexf(dataTX[:,1],dataTX[:,2:3]);

  DESCRIPTION="Exact solution for th0=cos(2*pi*y); kappa=0; v=[0,sin(2*pi*x)]";
  HOWTO="thexf(t,x) = 0.5+0.5*cos( 2*pi*(x[:,2]-sin(2*pi*x[:,1]).*t) )";
end

## MAKE DATA ##

noise = zeros(size(G_true));  @printf("not adding noise for now...\n");
y = G_true + noise;

## WRITE RESULTS TO FILE ##

#write
if (h5fl!="none")
  #first check to see if file exists and change names if it does
  h5fl0=h5fl; 
  i=0; 
  h5fl=replace(h5fl0,".h5",@sprintf("_%03d.h5",i)); 
  while (isfile(h5fl)) 
    @printf "file %s already exists\n" h5fl; 
    i+=1; 
    h5fl=replace(h5fl0,".h5",@sprintf("_%03d.h5",i)); 
  end 
  @printf("using file %s\n",h5fl);

  fid = h5open(h5fl,"w");
  isdefined(:knorm) && write(fid, "knorm", knorm );  #
  write(fid, "kappa" , kappa   );  #Save kappa
  write(fid, "t0"    , ad.t0   );  #Save t0
  write(fid, "dt"    , ad.dt   );  #Save dt
  write(fid, "tf"    , ad.tf   );  #Save tf
  write(fid, "th0k"  , thk2rl(th0k)    );  #Save initial condition
  isdefined(:vk_true) && write(fid, "vk_true", vk2rl(vk_true) );  #Save vector field used to generate measurement data
  isdefined(:v_true)  && write(fid, "v_true", v       );          #Save vector field used to generate measurement data
  isdefined(:DESCRIPTION)  && write(fid, "DESCRIPTION", DESCRIPTION   );  #Save description of purpose of this data file
  isdefined(:HOWTO)  && write(fid, "HOWTO", HOWTO   );                    #Save description of how to remake this data file
  write(fid, "dataTX", dataTX  );  #Save measurement points 
  write(fid, "G_true" , G_true );  #Save G(v) (data absent noise)
  write(fid, "y"      , y      );  #Save measurement data
  write(fid, "noise"  , noise  );  #Save measurement noise
  close(fid);
  @printf "Data written to HDF5 archive: %s\n" h5fl;
else
  @printf "Data not written to HDF5 archive (h5fl=\"%s\").\n" h5fl;
end

