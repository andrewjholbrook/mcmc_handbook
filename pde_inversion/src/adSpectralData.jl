#Data types for observations of the passive scalar for the advection-diffusion equation

#Data type for spectral data
mutable struct adSpectralData

  datatype::String
  dim::Int64        #Number of observations

  #tobs,kobs such that obs(x)=x[tobs,kobs]
  tobs::Array
  kobs::Array

  y::Array          #Data

  function adSpectralData()
    d = new()
    d.datatype = "spectral"
    return d;
  end
end
