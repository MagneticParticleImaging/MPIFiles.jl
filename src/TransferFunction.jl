export TransferFunction, sampleTF, tf_receive_chain

mutable struct TransferFunction
  freq::Vector{Float64}
  data::Matrix{ComplexF64}
  inductionFactor::Vector{Float64}

  function TransferFunction(freq_, datain::Array{T}, inductionFactor=ones(size(datain,2))) where {T<:Complex}
    freq = freq_[1]:(freq_[2]-freq_[1]):freq_[end]
    data=reshape(deepcopy(datain),size(datain,1), size(datain,2))
    return new(freq_, data, inductionFactor)
  end
end

function TransferFunction(freq_, ampdata, phasedata, args...)
  data = ampdata.*exp.(im.*phasedata)
  return TransferFunction(freq_, data, args...)
end

function TransferFunction(filename::String; kargs...)
    filenamebase, ext = splitext(filename)
    if  ext == ".h5"
      tf = load_tf(filename)
    else #if ext == "s1p" || ext == "s2p"
      tf = load_tf_fromVNA(filename; kargs...)
    end
    return tf
end

function getindex(tmf::TransferFunction, x::UnitRange, chan::Integer=1)
  a = tmf.data[x,chan]
  return a
end

#function getindex(tmf::TransferFunction, x::Real, chan::Integer=1)
#  a = tmf.interp[chan](x)
#  return a
#end

function getindex(tmf::TransferFunction, X::Union{Vector,AbstractRange}, chan::Integer=1)
  I = extrapolate(interpolate((tmf.freq,), tmf.data[:,chan], Gridded(Linear())), Interpolations.Flat())

  return [I(x) for x in X]
end

function getindex(tmf::TransferFunction, X::Union{Vector,AbstractRange}, chan::AbstractRange)
  out = zeros(ComplexF64, length(X), length(chan))
  for d=1:length(chan)
    out[:,d] = tmf[X,d]
  end
  return out
end

function load_tf(filename::String)
  tf = h5read(filename,"/transferFunction")
  freq = h5read(filename,"/frequencies")
  inductionFactor = h5read(filename,"/inductionFactor")
  return TransferFunction(freq,tf,inductionFactor)
end

function combine(tf1, tf2)
  freq = tf1.freq
  data = cat(tf1.data,tf2.data, dims=2)
  inductionFactor = cat(tf1.inductionFactor, tf2.inductionFactor, dims=1)
  return TransferFunction(freq, data, inductionFactor)
end

function save(filename::String, tf::TransferFunction)
  h5write(filename, "/transferFunction", tf.data)
  h5write(filename, "/frequencies", tf.freq)
  h5write(filename, "/inductionFactor", tf.inductionFactor)
  return nothing
end

function load_tf_fromVNA(filename::String;
    frequencyWeighting=false,
    R = 50.0, #Ω
    N = 10, #5# Turns
    A = 7.4894*10.0^-4) #m^2 #1.3e-3^2*pi;)

  file = open(filename)
  lines = readlines(file)
  apdata = Float64[]
  aϕdata = Float64[]
  freq = Float64[]
  for i=4:length(lines)
      tmp = split(strip(lines[i])," ")
      tmp=tmp[tmp.!=""]
      f = parse(Float64,strip(tmp[1]))
      if occursin(lines[3],"# kHz S MA R 50")
          ap = parse(Float64,strip(tmp[2]))
          aphi = parse(Float64,strip(tmp[3]))
          f=f*1000
       elseif occursin(lines[3],"# Hz S DB R 50")
            ap = 10^(parse(Float64,strip(tmp[2]))/20)
            aphi = parse(Float64,strip(tmp[3]))*pi/180
            f=f
       elseif occursin(lines[3],"# kHz S DB R 50")
            ap = 10^(parse(Float64,strip(tmp[2]))/20)
            aphi = parse(Float64,strip(tmp[3]))*pi/180
            f=f*1000
        elseif occursin(lines[3],"# MHz S DB R 50")
            ap = 10^(parse(Float64,strip(tmp[2]))/20)
            aphi = parse(Float64,strip(tmp[3]))*pi/180
            f=f*1000000
      elseif occursin(lines[3],"# kHz S RI R 50")
          tf_complex=parse(Float64,strip(tmp[2]))+im*parse(Float64,strip(tmp[3]))
          ap=abs.(tf_complex);
          aphi=angle(tf_complex)
          f=f*1000
      else
	      error("Wrong data Format! Please export in kHz domain S21 parameter with either Magnitude/Phase, DB/Phase or Real/Imaginary!")
      end
      push!(apdata, ap)
      push!(aϕdata, aphi)
      push!(freq, f)
  end
  close(file)
  if frequencyWeighting
  	apdata ./= (freq.*2*pi) # As TF is defined as u_ADC = u_coil *TF the derivative from magnetic moment is applied as the division of the TF by w
  end
  return TransferFunction(freq, apdata, aϕdata)
end


function sampleTF(tmf::TransferFunction, f::MPIFile)
  freq = rxFrequencies(f)
  numChan = rxNumChannels(f)
  numFreq = length(freq)
  return tmf[freq,1:numChan]
end
