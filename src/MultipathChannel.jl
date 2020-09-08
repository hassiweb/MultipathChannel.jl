module MultipathChannel

using FractionalDelayFilter
using LinearAlgebra
using DSP
using Statistics

export ChanPars, initChan, applyChan

"""
    MultipathChannel.ChanPars
Type of channel parameters containing:

- `numpaths::Int`: Number of paths
- `pathdelays::Vector{Float64}`: Delay of each path in second.  This parameter is set according to `channelmodel` when calling `MultipathChannel.initChan`.
- `averagepathgains::Vector{Float64}`: Average gain of each path in dB.  This parameter is set according to `channelmodel` when calling `MultipathChannel.initChan`.
- `fdfcoefs::Vector{Vector{Float64}}`: Coefficients of fractional delay filter (FDF).  This parameter is obtained by `pathdelays`.
- `nonzeroindices::Vector{Int}`: Non zero index of each path for FDF.  This parameter is obtained by `pathdelays`.
- `pathdopplers::Vector{Float64}`: Doppler frequency of each path in Hz.  Doppler frequency is given according to Jakes model with `maxdoppler`.  This parameter is initialized when calling `initChan`.
- `fadingamplitudes::Vector{Complex{Float64}}`: Instantaneous gain of each path given by Rayleigh fading.  This parameter is initialized when calling `initChan`.
- `initialphases::Vector{Float64}`: Initial phase offset of each path.  This parameter is initialized when calling `initChan`.
- `samplingrate::Float64`:  Sampling rate in Hz
"""
struct ChanPars
    numpaths::Int
    pathdelays::Vector{Float64} 
    averagepathgains::Vector{Float64}
    fdfcoefs::Vector{Vector{Float64}}
    nonzeroindices::Vector{Int}
    pathdopplers::Vector{Float64}
    fadingamplitudes::Vector{Complex{Float64}}
    initialphases::Vector{Float64}
    samplingrate::Float64
end

"""
    MultipathChannel.initChan(maxdoppler, samplingrate, channelmodel; rmsdelay, maxorder, gainnormalization=true)
Generate a `ChanPars` type channel parameters.

**Arguments**
- `maxdoppler`: Maximum Doppler frequency (Hz)
- `samplingrate`: Sampling rate (Hz)
- `channelmodel`: Channel model for path delays and average path gains: supporting TDL-C, and single Rayleigh path
- `rmsdelay`: Root mean square delay for TDL channel models in nano seconds
- `maxorder`: Maximum filter order for FDF
- `gainnormalization`: Indicator for normalizing path gains
"""
function initChan(maxdoppler, samplingrate, channelmodel; rmsdelay=300.0, maxorder=8, gainnormalization=true)
    if lowercase(channelmodel) == "tdlc"
        pathdelays = rmsdelay * samplingrate * [0.000; 0.2099; 0.2219; 0.2329; 0.2176; 0.6366; 0.6448; 0.6560; 0.6584; 0.7935; 0.8213; 0.9336; 1.2285; 1.3083; 2.1704; 2.7105; 4.2589; 4.6003; 5.4902; 5.6077; 6.3065; 6.6374; 7.0427; 8.6523]
        averagepathgains = [-4.4; -1.2; -3.5; -5.2; -2.5; 0; -2.2; -3.9; -7.4; -7.1; -10.7; -11.1; -5.1; -6.8; -8.7; -13.2; -13.9; -13.9; -15.8; -17.1; -16; -15.7; -21.6; -22.8]
    elseif lowercase(channelmodel) == "rayleigh"
        pathdelays = [0.00]
        averagepathgains = [0.00]  # in dB
    end

    if gainnormalization
        averagepathgains = 10*log10.( @.(10^(averagepathgains/10)) / sum(@.(10^(averagepathgains/10))) )
    end

    numpaths = length(averagepathgains)

    # Initialize Doppler frequencies, initial phase offsets, and fading amplitudes
    aoa = 2π*rand(numpaths)  # angle of arrival
    pathdopplers = maxdoppler * cos.(aoa)
    initialphases = 2π*rand(numpaths)
    fadingamplitudes = (randn(numpaths)+1im*randn(numpaths)) / √2  # Rayleigh fading

    # Initialize fractional delay filter
    fdfcoefs = Vector{Vector{Float64}}(undef, numpaths)
    nonzeroindices = Vector{Int}(undef, numpaths)
    filterorders = Vector{Int}(undef, numpaths)
    for n = 1:numpaths
        filterorders[n] = filtord(pathdelays[n], maxorder=maxorder)
        (fdfcoefs[n], nonzeroindices[n]) = getfdfcoef(filterorders[n], pathdelays[n])
    end

    return ChanPars(
        numpaths, 
        pathdelays, 
        averagepathgains, 
        fdfcoefs, 
        nonzeroindices, 
        pathdopplers, 
        fadingamplitudes, 
        initialphases, 
        samplingrate
    )
end

"""
    MultipathChannel.applyChan(in, chan::ChanPars; initialtime)
Apply multipath channel based on `chan`.

**Arguments**
- `in`: Input signal
- `chan::ChanPars`: Channel parameters
- `initialtime`: Time for generating Doppler shift
"""
function applyChan(in, chan::ChanPars; initialtime=0.0)
    delayedwave = zeros(ComplexF64, length(in)+length(chan.fdfcoefs[end])-2+chan.nonzeroindices[end], 1)
    rayleighwave = zeros(ComplexF64, length(in)+length(chan.fdfcoefs[end])-2+chan.nonzeroindices[end], 1)
    t = initialtime .+ (0:length(delayedwave[:,1])-1)/chan.samplingrate
    for n = 1:chan.numpaths
        fdfilter!(view(delayedwave, 1:length(in)+length(chan.fdfcoefs[n])-2+chan.nonzeroindices[n]), in, chan.fdfcoefs[n], chan.nonzeroindices[n])
        @. rayleighwave += sqrt(10^(chan.averagepathgains[n]/10)) * chan.fadingamplitudes[n] * exp(1im*2π*(chan.pathdopplers[n]*t-chan.initialphases[n])) * delayedwave
    end
    return rayleighwave
end

end