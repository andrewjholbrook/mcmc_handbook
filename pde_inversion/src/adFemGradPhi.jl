#adFemGradPhi() computes the gradient of phi by computing 
# rho, the solution to the adjoint problem, then
# D_{e_k} \Phi = \int_0^T \int_D \rho*(e_k \cdot \nabla \theta)
#for each e_k
#See Section 3.3 of https://arxiv.org/pdf/1808.01084.pdf

# inputs:
#  ad        structure describing advection-diffusion problem
#  kdisc     structure describing spectral discretization
#  mesh      structure describing fem mesh
#  theta     solution to the advection-diffusion equation, array dim=(nt,nUnk)
#  rho       solution to the adjoint equation, array dim=(nt,nUnk)
#   
# outputs:
#  gradPhi   gradient of \Phi wrt components of v
#
function adFemGradPhi(ad::adProb, kdisc::specDiscrete2D,mesh::mesh2dPbc,theta::Array,dPhidG::Array,data)

#  error("adFemGradPhi() does not work yet. need to figure out how to iterate through directions vDir for each basis.\n");

  #Assemble adjoint forcing
  if data.datatype == "point"
    ##DELTA FUNCTION VERSION
    ##Create delta functions at (t_j,x_j)
    ##  F[i] = \sum_j \Sigma^{-1} [y_j - \theta(t_j,x_j)] \phi_j
    ##  iso does interpolation (of phi) in x
    ##  tmat does interpolation (of theta) in t
    #F=zeros(ad.nt,mesh.nUnk);
    ##loop through data points and nodes and add to F
    #for i=1:dataDim;
    #  for j=1:length(mesh.eConn[idx[i],:]); 
    #    F[:,mesh.nodesToUnk[mesh.eConn[idx[i],j]]]+=(-dPhidG[i]./ad.dt)*data.iso[i,j]*data.tmat[i,:]; 
    #  end
    #end
    F = sum([-dPhidG[i]*data.L2Kern[i] for i=1:length(data.L2Kern)]);
  #elseif data.datatype == "spectral"
  #
  else
    error("datatype $(data.datatype) unrecognized or unsupported.");

  end

  #Solve adjoint problem
  rho  = adFemAdjoint(ad, mesh, F);

  #time quadrature rule (trapezoidal)
  wt = trapRule(ad.t);

  #spatial quadrature rule
  rule  = 7;
  r,s,w = twodQuadratureRule(rule); 
  o     = ones(rule);

  nElDOF = size(mesh.eConn,2);

  xg  = Array{Float64}(undef,rule,2);
  wg  = Array{Float64}(undef,rule);
  phi = Array{Float64}(undef,rule,nElDOF);
  p_x = Array{Float64}(undef,rule,nElDOF);
  p_y = Array{Float64}(undef,rule,nElDOF);

  vDirg  = Array{Float64}(undef,rule,2);

  gradPhi = zeros(Complex{Float64},kdisc.nk); 


  for k=1:mesh.nElements
   
    nLocal    = mesh.eConn[k,:][:];
    xLocal    = mesh.x[nLocal,:];

    rhoLoc =   rho[:,mesh.nodesToUnk[nLocal]]; #rho at local nodes
    thLoc  = theta[:,mesh.nodesToUnk[nLocal]]; #theta at local nodes

    xg,wg,phi,p_x,p_y = twodShape( xLocal, r, s, w );
    vbg = vfbasis(kdisc,xg); #vector field basis functions at xg
  
    for l=1:kdisc.nk
      vDirg = vfXY(vbg,float((1:kdisc.nk).==l); forceReal=false);     # vector field evaluated at quadrature points

      # (A)_{ij} = <\phi_i,v_x d_x \phi_j> + <\phi_i,v_y d_y \phi_j>
      #ALoc   = twodBilinear( vDirg[:,1], phi, p_x, wg );
      #ALoc  += twodBilinear( vDirg[:,2], phi, p_y, wg );
      ALoc   = twodBilinear( vDirg[:,1], p_x, phi, wg );
      ALoc  += twodBilinear( vDirg[:,2], p_y, phi, wg );
      
      #integrate in time (for loop should be faster than outer product and diag())
      for i=1:ad.nt
        gradPhi[l] += wt[i] * dot(rhoLoc[i,:], ALoc * thLoc[i,:]);
      end
    end
  end

  return gradPhi;
end
