#Data type for scalar variance data
mutable struct adScalVarData

  datatype::String
  dim::Int64        #Number of observations

  y::Array          #Data

  t::Array          #Time of observations
  tmat::Array       #[dim,nt] interpolation matrix

  function adScalVarData()
    d = new()
    d.datatype = "scalvar"
    return d;
  end
  function adScalVarData(y)
    d = adScalVarData()
    d.dim    = length(y);
    #d.t      = t;
    d.y      = y;
    #d.tmat   = timeInterp(ad.t,dataTX[:,1]);
    return d;
  end
  #tobs: time of observations
  #tarray: timesteps of solver
  function adScalVarData(y,tobs,tarray)
    d = adScalVarData(y)
    d.t      = tobs;
    d.tmat   = timeInterp(tarray,tobs);
    return d;
  end
end
