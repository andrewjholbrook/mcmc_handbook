#kraichnanEnergy() computes the energy spectrum given in chen1998simulations equation (28)
# This version computes energies for k=1:kMax
# Inputs:
#  kMax       maximum k
#  zetaInf    \zeta_2^\infty; default value from Table I
#  nFields    number of fields (N); default value from 8192^2 simulation
#
# Outputs:
#  E          energies for k=1:kMax
#
function kraichnanEnergy(kMax::Int64; zetaInf=0.5, nFields=20)
  
  E = zeros(kMax); 

  A=1; #scaling factor
  
  #mentioned in between equations (28) & (29)
  xi = 2-zetaInf;
  
  #equation (28)
  for k=1:kMax
    for j=0:nFields
      kj = sqrt(2)^j;
      E[k] += A*( (k/kj)^4 * exp(-1.5*(k/kj)^2) * (kj^(-xi)) );
    end
  end

  return E;
end

#kraichnanEnergy() computes the energy spectrum given in chen1998simulations equation (28)
# This version computes energies for a list of k's
# Inputs:
#  kList      list of k's
#  zetaInf    \zeta_2^\infty; default value from Table I
#  nFields    number of fields (N); default value from 8192^2 simulation
#  A          scaling factor fixing total energy; default is 1
#
# Outputs:
#  E          energies
#
function kraichnanEnergy(kList::Array; zetaInf=0.5, nFields=20, A=1)
  
  nk = length(kList);

  E = zeros(nk); 

  #mentioned in between equations (28) & (29)
  xi = 2-zetaInf;
  
  #equation (28)
  for i=1:nk
    for j=0:nFields
      k = kList[i];
      kj = sqrt(2)^j;
      E[i] += A*( (k/kj)^4 * exp(-1.5*(k/kj)^2) * (kj^(-xi)) );
    end
  end

  @printf("scaling kraichnan prior so that radial integral matches chen1998simulations figure 6/7\n");
  #E[kList.!=0.0] = E[kList.!=0.0] ./ (2*pi*kList[kList.!=0.0]);
  for i=1:nk
    if kList[i] != 0.0
      E[i] = E[i] ./ (2*pi*kList[i]);
    end
  end

  return E;
end
#kraichnanEnergy() with kList as Range rather than Array
function kraichnanEnergy(kList::AbstractRange; zetaInf=0.5, nFields=20)
  return kraichnanEnergy(collect(kList); zetaInf=zetaInf, nFields=nFields);
end
