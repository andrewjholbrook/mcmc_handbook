#Computes the gradient of phi by computing 
# rho, the solution to the adjoint problem, then
# D_{e_k} \Phi = \int_0^T \int_D \rho*(e_k \cdot \nabla \theta)
#for each e_k
#See Section 3.3 of https://arxiv.org/pdf/1808.01084.pdf
#
# inputs:
#  ad        structure describing advection-diffusion problem
#  kdisc     structure describing spectral discretization
#  theta     solution to advection-diffusion problem
#   
# outputs:
#  gradPhi   gradient of \Phi wrt components of v
#
function adSpectralGradPhi(ad::adProb,kdisc::specDiscrete2D,theta::Array,dPhidG::Array,data)
  
  #Assemble adjoint forcing
  if data.datatype == "point"
    #dPhidG = dPhi(Gj);
 
    #spectral representation of delta function
    #dk = deltak(kdisc,data.dataTX[:,2:3]);
    #F = -data.tmat'*diagm(0 => dPhidG)*dk/ad.dt;
    F = sum([-dPhidG[i]*data.L2Kern[i] for i=1:length(data.L2Kern)]);

  elseif data.datatype == "spectral"
    #dPhidG = dPhi(Gj);
    F = zeros(ad.nt,kdisc.nk);
    F[data.tobs,data.kobs] = reshape(-dPhidG,length(data.tobs),length(data.kobs));

  elseif data.datatype == "scalvar"
    #dPhidG = dPhi(Gj);
    F = copy(theta);
    F[:,kdisc.kn == 0.0] = 0.0; #mean 0
    for i=1:size(theta,1)
      F[i,:] *= -2.0*dPhidG[i];
    end

  #elseif data.datatype=="test"
  #  #dPhidG = -Q0invcov()*(y - Gj)[:];
  #  dPhidG = dPhi(Gj,y);
  #  F = zeros(Complex{Float64},ad.nt,kdisc.nk);
  #  for j=1:dataDim
  #    Kj = Kjf(j);
  #    F += -dPhidG[j] * Kj[ad.nt:-1:1,:];
  #  end
  else
    error("adSpectralGradPhi: datatype $(data.datatype) unrecognized or unsupported.");

  end


  #Solve adjoint problem
  rho  = adSpectralAdjoint(ad, kdisc, F);

  #time quadrature rule (trapezoidal)
  wt = trapRule(ad.t);

  #precompute the dot product of all kperp times all k
  kpTimesK = kdisc.kperp*transpose(kdisc.k);
  #compute the time integral of all rho components against all theta components
  intRhoTh = transpose(rho)*diagm(0 => wt)*theta;

  if kdisc.complex
    gradPhi = zeros(Complex{Float64},kdisc.nk); 
  
  # \int_{0}^{T} \int_D \rho \left( \tilde{\vfield} \cdot \nabla \pdesol \right) = \sum_l \sum_j 2 \pi i \tilde{v}_l \left( \kbf_l^{\perp} \cdot \kbf_j \right) \int_{0}^{T} \rho_i(t) \pdesol_j(t) \quad\text{where}\quad i \ni \kbf_i = -\kbf_l-\kbf_j\\
  #  for l=1:kdisc.nk
  #    for j=1:kdisc.nk
  #      i=kdisc.kdiff[rev(j),l];
  ##      if (kdisc.k[l,:]+kdisc.k[j,:]+kdisc.k[i,:]) != zeros(2)
  ##        @printf("we didn't find the right i for l=%d, j=%d\n",l,j);
  ##      end
  #      if i!=0
  #        gradPhi[l] += 2*pi*im*dot(kdisc.kperp[l,:],kdisc.k[j,:])*dot(wt,rho[:,i].*theta[:,j]);
  #      end
  #    end
  #  end
    for l=1:kdisc.nk
      for j=1:kdisc.nk
        i=kdisc.kdiff[rev(j),l];
  #      if (kdisc.k[l,:]+kdisc.k[j,:]+kdisc.k[i,:]) != zeros(2)
  #        @printf("we didn't find the right i for l=%d, j=%d\n",l,j);
  #      end
        if i!=0
          #gradPhi[l] += 2*pi*im*dot(kdisc.kperp[l,:],kdisc.k[j,:])*dot(wt,rho[:,i].*theta[:,j]);
          gradPhi[l] += 2*pi*im*kpTimesK[l,j]*intRhoTh[i,j];
        end
      end
    end
  
  else #real
    error("The real version of the generic adSpectralGradPhi has not been implemented. Please use the templated version.");
  end #end real

  return gradPhi;
end

function adSpectralGradPhi(ad::adProb,kdisc::specDiscrete2D,Ip, Tp, In, Tn, T0, theta::Array,dPhidG::Array,data)
  
  #Assemble adjoint forcing
  if data.datatype == "point"
    #dPhidG = dPhi(Gj);
 
    #spectral representation of delta function
    #dk = deltak(kdisc,data.dataTX[:,2:3]);
    #F = -data.tmat'*diagm(0 => dPhidG)*dk/ad.dt;
    F = sum([-dPhidG[i]*data.L2Kern[i] for i=1:length(data.L2Kern)]);

  elseif data.datatype == "spectral"
    #dPhidG = dPhi(Gj);
    F = zeros(ad.nt,kdisc.nk);
    F[data.tobs,data.kobs] = reshape(-dPhidG,length(data.tobs),length(data.kobs));

  elseif data.datatype == "scalvar"
    #dPhidG = dPhi(Gj);
    F = copy(theta);
    F[:,kdisc.kn == 0.0] = 0.0; #mean 0
    for i=1:size(theta,1)
      F[i,:] *= -2.0*dPhidG[i];
    end

  #elseif data.datatype=="test"
  #  #dPhidG = -Q0invcov()*(y - Gj)[:];
  #  dPhidG = dPhi(Gj,y);
  #  F = zeros(Complex{Float64},ad.nt,kdisc.nk);
  #  for j=1:dataDim
  #    Kj = Kjf(j);
  #    F += -dPhidG[j] * Kj[ad.nt:-1:1,:];
  #  end
  else
    error("adSpectralGradPhi: datatype $(data.datatype) unrecognized or unsupported.");

  end

  return adSpectralGradPhi(ad, kdisc, Ip, Tp, In, Tn, T0, theta, F);
end

#in this version the forcing is provided explicitly
function adSpectralGradPhi(ad::adProb, Ip, Tp, In, Tn, T0, theta::Array{Float64,2},F::Array{Float64,2})

  #Solve adjoint problem
  rho  = adSpectralAdjoint(ad, Ip, Tp, In, Tn, T0, F);

  #time quadrature rule (trapezoidal)
  wt = trapRule(ad.t);

  ##precompute the dot product of all kperp times all k
  #kpTimesK = kdisc.kperp*transpose(kdisc.k);
  #compute the time integral of all rho components against all theta components
  intRhoTh = transpose(rho)*diagm(0 => wt)*theta;

  #if kdisc.complex
  #  error("The complex version of the templated adSpectralGradPhi has not been implemented. Please use the generic version.");
  #else #real
  #end #end real

  #now that we have rho, we need to compute \int \rho ek \cdot \theta
  #this is the same as -rho'*A(ek)*theta where A(ek) is the stiffness matrix associated with the advection term of the AD equation. Ip and In tell us exactly how those terms map to components ek of v, so we iterate across them and add to the relevant terms in the gradient.
  gradPhi = zeros(length(ad.v)); 
  Mp = -Tp .* intRhoTh;
  Mn = -Tn .* intRhoTh;
  
  for i=1:size(Ip,1)
    gradPhi[ Ip[i,2] ] += Mp[Ip[i,1]];
  end 
  for i=1:size(In,1)
    gradPhi[ In[i,2] ] += Mn[In[i,1]];
  end 

  return gradPhi;
end
