#Data types for observations of the passive scalar for the advection-diffusion equation

#Data type for point data
mutable struct adPointData

  datatype::String
  dim::Int64        #Number of observations
  nPts::Int64       #Number of observation points

  y::Array          #Data

  dataTX::Array     #[dim,3] array of (t,x,y) observation locations 
  tmat::Array       #[dim,nt] interpolation matrix

#  basisXY            #Spatial basis functions evaluated at x,y observation points

  #These are only used by FEM
  idx::Array
  iso::Array

  L2Kern            #Kernel such that integrating against it in space and time produces the observations

  function adPointData()
    d = new()
    d.datatype = "point"
    return d;
  end
  function adPointData(dataTX,y)
    d = adPointData()
    d.nPts   = size(dataTX,1);
    d.dim    = length(y);
    d.dataTX = dataTX;
    d.y      = y;
    return d;
  end

  #here we know that they're using a spectral discretization
  function adPointData(ad::adProb,kdisc::specDiscrete2D,dataTX,y; computeL2Kern=true)
    d = adPointData(dataTX,y)

    d.tmat   = timeInterp(ad.t,dataTX[:,1]);

    #d.basisXY= e2piikx(kdisc.k,data.dataTX[:,2:3]);
    dk = deltak(kdisc,d.dataTX[:,2:3]);
    if computeL2Kern
      #d.L2Kern = [ d.tmat[i,:]*transpose(dk[i,:]/ad.dt) for i=1:size(dataTX,1) ];
      d.L2Kern = [ sparse(d.tmat[i,:]*transpose(dk[i,:]/ad.dt)) for i=1:size(dataTX,1) ];
    end

    return d;
  end

  #here we know that they're using a FEM discretization
  function adPointData(ad::adProb,mesh::mesh2dPbc,dataTX,y)
    d = adPointData(dataTX,y)

    d.tmat   = timeInterp(ad.t,dataTX[:,1]);

    #d.basisXY= e2piikx(kdisc.k,data.dataTX[:,2:3]);
    d.idx,d.iso = twodMeshSearch(mesh.eConn,mesh.x,dataTX[:,2:3]);

    d.L2Kern=[ spzeros(ad.nt,mesh.nUnk) for i=1:d.dim ];
    #for i=1:length(d.L2Kern)
    #  d.L2Kern[i][:,mesh.nodesToUnk[mesh.eConn[d.idx[i],:]]]+=d.tmat[i,:]*transpose(d.iso[i,:])/ad.dt;
    #end
    M = twodMassMatrix(mesh);
    for i=1:d.dim;
      phiX=zeros(mesh.nUnk);
      phiX[mesh.nodesToUnk[mesh.eConn[d.idx[i],:]]]+=d.iso[i,:];
      delX = M \ phiX;
      d.L2Kern[i]=d.tmat[i,:]*transpose(delX)/ad.dt;
      #o3[i]=dot(wt,twodIntegrate(mesh,transpose(L2Kern1).*transpose(th)))
    end

    return d;
  end
end

