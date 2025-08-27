#h5WriteStruct(): Save a struct to an hdf5 file
function h5WriteStruct(filename::String,data;verbose=1)
  f = h5open(filename,"w");
  write(f,"typeof",String(Symbol(typeof(data))));
  for field in fieldnames(typeof(data))
    write(f,String(field),getproperty(data,field));
    (verbose>1) && println("+ $(String(field))");
  end
  close(f);
  (verbose>0) && println("Wrote: $(filename)");
end

#h5ReadStruct(): Read a struct from an hdf5 file
#The struct type is read from the "typeof" field in the file
function h5ReadStruct(filename;verbose=1)
  f = h5open(filename,"r");
  #create an object of type read from the "typeof" field
  #(this may be far from the best way to do this)
  to = read(f,"typeof");
  s = eval(Symbol(to))();
  #iterate over fields of the type and read everything in
  for field in fieldnames(typeof(s))
    setproperty!(s,field,read(f,String(field)));
    (verbose>1) && println("+ $(String(field))");
  end
  close(f);
  (verbose>0) && println("Read from: $(filename)");
  return s;
end

#HDF5.write(): Extend to SparseCSC
function HDF5.write(f::HDF5.File,field::String,m::SparseMatrixCSC)
  (II,JJ,ZZ) = findnz(m);
  write(f,"$(field)/typeof","SparseMatrixCSC");
  #write(f,"$(field)/typeof",typeof(m));
  write(f,"$(field)/II",II);
  write(f,"$(field)/JJ",JJ);
  write(f,"$(field)/ZZ",ZZ);
end

#convert(): Extend to SparseCSC
function Base.convert(::Type{SparseMatrixCSC},d::AbstractDict)
  return sparse(d["II"],d["JJ"],d["ZZ"]);
end
