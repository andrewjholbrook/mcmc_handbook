#using Pkg;
#pkg"activate ../..";

using CUDA
using Flux
using Flux: @epochs, train!
using BSON: @save, @load
using ProgressMeter
using HDF5
using StatsBase
using Statistics
using Printf

#Parse command line arguments
using ArgParse
apSettings = ArgParseSettings();
@add_arg_table apSettings  begin
  "--modeltype"
    help = "model type"
    default = "original"
  "--epochs"
    help = "number of epochs"
    arg_type = Int
    default = 100
  "--trainpct"
    help = "% of samples to use for training"
    arg_type = Float64
    default = 0.8
end
args = parse_args(apSettings,as_symbols=true);
#display(args)

useGPU = true;


nEpochs    = args[:epochs];
modelType  = args[:modeltype];
trainRatio = args[:trainpct];

println("Model type: $(modelType)");
println("Epochs: $(nEpochs)");
println("Training with fraction of samples: $(trainRatio)");


ADR_ROOT = "/projects/SIAllocation"
dataFile      = ADR_ROOT*"/lolmc/juqex2_traindata.h5";
#modelFileBase = ADR_ROOT*"/lolmc/trained_models/model";
#metaFileBase  = ADR_ROOT*"/lolmc/trained_models/meta";
modelFileDir  = ADR_ROOT*"/lolmc/trained_models/"*modelType;
modelFileBase = modelFileDir*"/model";
metaFileBase  = modelFileDir*"/meta";

if ( ! isdir(modelFileDir) ) 
  println("Output directory $(modelFileDir) does not exist. Creating...");
  mkdir(modelFileDir);
end

#get model output file name
mkfilename(cnt;modelFileBase=modelFileBase) = @sprintf("%s_%03d.bson",modelFileBase,cnt);
mkmetaname(cnt;metaFileBase =metaFileBase ) = @sprintf("%s_%03d.h5"  ,metaFileBase ,cnt);
cnt=0;
while (isfile(mkfilename(cnt)))
  global cnt;
  cnt+=1;
end
modelOutFile = mkfilename(cnt);
metaOutFile  = mkmetaname(cnt);
#println("Model file: $(modelOutFile)");

#get model input file name
modelInFile = "";
metaInFile  = "";
if cnt > 0
  modelInFile = mkfilename(cnt-1);
  metaInFile  = mkmetaname(cnt-1);
end

#read data
f = h5open(dataFile);
samples = read(f,"samples");
gradPots = read(f,"gradPots");
close(f);

#drop first two columns
samples  = samples[:,3:end];
gradPots = gradPots[:,3:end];

nsamp, ndim = size(samples);

#load or build model & metadata
if modelInFile != ""
  println("Reading model from: $(modelInFile)");
  @load modelInFile model

  println("Reading metadata from: $(metaInFile)");
  f = h5open(metaInFile,"r");
  trainIdx = read(f,"trainIdx");
  ntrain   = length(trainIdx);
  testIdx  = read(f,"testIdx" );
  close(f);
else
  #model = Chain(
  #  BatchNorm(ndim),
  #  Dense( ndim,  588, softplus),
  #  Dense(  588,  980, relu),
  #  Dense(  980, 1960, softplus),
  #  Dense( 1960,  980, selu),
  #  Dense(  980,  588, softplus),
  #  Dense(  588, ndim )
  #  );
  if modelType == "original"
    model = Chain(
      BatchNorm(ndim),
      Dense( ndim,   512, softplus),
      Dense(  512,  1024, relu),
      Dense( 1024,  2048, softplus),
      Dense( 2048,  1024, selu),
      Dense( 1024,   512, softplus),
      Dense(  512,  ndim )
    );
  elseif modelType == "large"
    model = Chain(
      BatchNorm(ndim),
      Dense( ndim,   512, softplus),
      Dense(  512,  1024, relu),
      Dense( 1024,  2048, softplus),
      Dense( 2048,  4096, selu),
      Dense( 4096,  2048, softplus),
      Dense( 2048,  1024, selu),
      Dense( 1024,   512, softplus),
      Dense(  512,  ndim )
    );
  elseif modelType == "small"
    model = Chain(
      BatchNorm(ndim),
      Dense( ndim,   512, softplus),
      Dense(  512,  1024, relu),
      Dense( 1024,   512, softplus),
      Dense(  512,  ndim )
    );
  elseif modelType == "relusmall"
    model = Chain(
      BatchNorm(ndim),
      Dense( ndim,   512, relu),
      Dense(  512,  1024, relu),
      Dense( 1024,   512, relu),
      Dense(  512,  ndim )
    );
  elseif modelType == "selusmall"
    model = Chain(
      BatchNorm(ndim),
      Dense( ndim,   512, selu),
      Dense(  512,  1024, selu),
      Dense( 1024,   512, selu),
      Dense(  512,  ndim )
    );
  else
    println("ERROR: Unrecognized model type: $(modelType)");
    exit(1);
  end

  #indices of training and test data
  ntrain = round(Int,trainRatio*nsamp);
  trainIdx = sample(1:nsamp,ntrain,replace=false);
  testIdx  = setdiff(1:nsamp,trainIdx);
end

#get training data
datax = samples[trainIdx,:];
datay = gradPots[trainIdx,:];

# #normalization function
# #(don't use Flux's normalise() because we might need to reconstruct)
# function normalize( data )
#   dmean   = mean(data,dims=1)[:];
#   dstd    = std(data,dims=1)[:];
#   dnormal = ( data .- dmean' )./(dstd');
#   return dnormal, dmean, dstd;
# end
# 
# #reconstruct dataset from normalized version, mean, and std dev
# function reconstruct( dnormal, dmean, dstd )
#   return ( dmean' .+ (dnormal .* dstd') );
# end
# 
# #normalize inputs and outputs
# datax, xmean, xstd = normalize(datax);
# #datay, ymean, ystd = normalize(datay);

#transpose
datax = Matrix(datax');
datay = Matrix(datay');

#normalize test inputs and outputs
testx = samples[testIdx,:];
testy = gradPots[testIdx,:];
#testx, _, _ = normalize(testx);
#testy, _, _ = normalize(testy);
testx = Matrix(testx');
testy = Matrix(testy');

#loss function
#lossFcn(x, y) = Flux.Losses.mae(model(x), y);
#lossFcn(x,y=datay) = mean(abs.(model(x) .- y)[:]);
lossFcn(x,y) = Flux.Losses.mse(model(x),y);


# GPU
if useGPU
  #datax = cu(datax);
  #datay = cu(datay);
  #testx = cu(testx);
  #testy = cu(testy);
  #model = fmap(cu, model);
  datax = datax |> gpu;
  datay = datay |> gpu;
  testx = testx |> gpu;
  testy = testy |> gpu;
  model = model |> gpu;
  #lossFcn(x) = mean(abs.(model(x) .- datay)[:]);
  #lossFcn(x,y=datay) = mean(abs.(model(x) .- y)[:]);
  #lossFcn(x,y=datay) = Flux.Losses.mse(x,y);

#  loader = Flux.DataLoader((datax, datay) |> gpu, batchsize=64, shuffle=true);
#else
end

#optimizer
optim = Flux.setup(Flux.Adam(), model); #(Flux.Adam(0.001, (0.9,0.999), 1.0e-8), model)

#data loader
loader = Flux.DataLoader((datax, datay), batchsize=64, shuffle=true);

lossTrainStart = lossFcn(datax,datay);
lossTestStart  = lossFcn(testx,testy);
println("Loss on training data is: $( lossTrainStart )");
println("Loss on test data is:     $( lossTestStart  )");

#train
println("Starting training.")
lossesTrain = [];
lossesTest  = [];
@showprogress for epoch in 1:nEpochs #1_000
  #Flux.train!(model, loader, optim) do m, x, y
  #  y_hat = m(x);
  #  Flux.Losses.mse(y_hat,y);
  #end
  for (x, y) in loader
    loss, grads = Flux.withgradient(model) do m
      # Evaluate model and loss inside gradient context:
      y_hat = m(x);
      Flux.Losses.mse(y_hat,y)
    end
    Flux.update!(optim, model, grads[1]);
    push!(lossesTrain, loss);  # logging, outside gradient context

    y_hat = model(testx);
    lossTest = Flux.Losses.mse(y_hat,testy)
    push!(lossesTest, lossTest);  # logging, outside gradient context
  end
end
println("Done training.");

lossesTrainEpoch = mean.(Iterators.partition(lossesTrain, length(loader)));
lossesTestEpoch  = mean.(Iterators.partition(lossesTest , length(loader)));

@printf("%8s %12s %12s\n","Epoch","Loss (Train)", "Loss (Test)");
for ep=1:10:nEpochs
  @printf("%8d %12.4f %12.4f\n",ep,lossesTrainEpoch[ep],lossesTestEpoch[ep]);
end

minTrainEpoch = argmin(lossesTrainEpoch);
minTestEpoch  = argmin(lossesTestEpoch);
println("Minimum training loss was $(minimum(lossesTrainEpoch)) at epoch $(minTrainEpoch)");
println("Minimum testing  loss was $(minimum(lossesTestEpoch )) at epoch $(minTestEpoch )");

lossTrainEnd = lossFcn(datax,datay);
lossTestEnd  = lossFcn(testx,testy);
println("Loss on training data is: $( lossTrainEnd )");
println("Loss on test data is:     $( lossTestEnd  )");

# #test on cpu to make sure everything works before saving
# println("Moving back to CPU:");
# datax = datax |> cpu;
# datay = datay |> cpu;
# testx = testx |> cpu;
# testy = testy |> cpu;
# model = model |> cpu;
# lossFcn(x,y=datay) = mean(abs.(model(x) .- y)[:]);
# println("Loss on training data is: $( lossFcn(datax,datay) )");
# println("Loss on test data is:     $( lossFcn(testx,testy) )");

#save model
model = model |> cpu;
@save modelOutFile model
println("Wrote: $(modelOutFile)");

#save metadata
f = h5open(metaOutFile,"w");
write(f,"dataFile",dataFile);
write(f,"trainRatio" ,trainRatio );
write(f,"trainIdx",trainIdx);
write(f,"testIdx" ,testIdx );
write(f,"nEpochs" ,nEpochs);
write(f,"lossTrainStart" ,lossTrainStart );
write(f,"lossTestStart"  ,lossTestStart  );
write(f,"lossTrainEnd"   ,lossTrainEnd   );
write(f,"lossTestEnd"    ,lossTestEnd    );
write(f,"lossesTrainBatch",Float64.(lossesTrain));
write(f,"lossesTestBatch" ,Float64.(lossesTest ));
write(f,"lossesTrainEpoch",lossesTrainEpoch);
write(f,"lossesTestEpoch" ,lossesTestEpoch );
if modelInFile != ""
  write(f,"modelInFile" ,modelInFile );
  write(f,"metaInFile" ,metaInFile );
end
write(f,"modelOutFile" ,modelOutFile );
close(f);
println("Wrote: $(metaOutFile)");

# using Plots
# plot(losses; xaxis=(:log10, "iteration"), yaxis="Loss", label="per batch");
# n = length(loader);
# plot!(n:n:length(losses), mean.(Iterators.partition(losses, n)),label="epoch mean", dpi=200);
# savefig(modelFileDir*"/loss.png");
# println("Wrote: "*modelFileDir*"/loss.png");
