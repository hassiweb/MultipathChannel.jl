module MultipathChannel

using FractionalDelayFilter

export ChanPars, initChan, applyChan

struct ChanPars
    numpaths::Int
    pathdelays::Vector{Float64}
    pathdopplers::Vector{Float64}
    pathgains::Vector{Float64}
    pathoffset::Vector{Float64}
    pathrayleighamp::Vector{Float64}
    fdfcoefs::Vector{Vector{Float64}}
    nonzeroindices::Vector{Int}
    samplingrate::Float64
end

const j = 1im

function initChan(maxdoppler, samplingrate, channelmodel; rmsdelay=300, maxorder=8)
    if lowercase(channelmodel) == "tdlc"
        pathdelays = rmsdelay * samplingrate * [0.000; 0.2099; 0.2219; 0.2329; 0.2176; 0.6366; 0.6448; 0.6560; 0.6584; 0.7935; 0.8213; 0.9336; 1.2285; 1.3083; 2.1704; 2.7105; 4.2589; 4.6003; 5.4902; 5.6077; 6.3065; 6.6374; 7.0427; 8.6523]
        pathgains = [-4.4; -1.2; -3.5; -5.2; -2.5; 0; -2.2; -3.9; -7.4; -7.1; -10.7; -11.1; -5.1; -6.8; -8.7; -13.2; -13.9; -13.9; -15.8; -17.1; -16; -15.7; -21.6; -22.8]
    elseif lowercase(channelmodel) == "rayleigh"
        pathdelays = [0.000; 0.000; 0.000; 0.000; 0.000]
        pathgains = [0.00; 0.00; 0.00; 0.00; 0.00]
    end

    numpaths = length(pathgains)
    pathdopplers = Vector{Float64}(undef, numpaths)
    pathoffsets = Vector{Float64}(undef, numpaths)
    pathrayleighamp = Vector{Float64}(undef, numpaths)
    for n = 1:numpaths
        aoa = 2π*rand()  # angle of arrival
        pathdopplers[n] = maxdoppler * cos(aoa)
        pathoffsets[n] = 2π*rand()
        pathrayleighamp[n] = √((randn()^2+randn()^2)/length(pathgains))
    end

    # initialize fractional delay filter
    fdfcoefs = Vector{Vector{Float64}}(undef, numpaths)
    nonzeroindices = Vector{Int}(undef, numpaths)
    filterorders = Vector{Int}(undef, numpaths)
    for n = 1:numpaths
        filterorders[n] = filtord(pathdelays[n], maxorder=maxorder)
        (fdfcoefs[n], nonzeroindices[n]) = getfdfcoef(filterorders[n], pathdelays[n])
    end

    return ChanPars(numpaths, pathdelays, pathdopplers, pathgains, pathoffsets, pathrayleighamp, fdfcoefs, nonzeroindices, samplingrate)
end

function applyChan(in, chan::ChanPars; initialtime=0.0)
    delayedwave = zeros(ComplexF64, length(in)+length(chan.fdfcoefs[end])-2+chan.nonzeroindices[end], 1)
    rayleighwave = zeros(ComplexF64, length(in)+length(chan.fdfcoefs[end])-2+chan.nonzeroindices[end], 1)
    t = initialtime .+ (0:length(delayedwave[:,1])-1)/chan.samplingrate
    for n = 1:chan.numpaths
        fdfilter!(view(delayedwave, 1:length(in)+length(chan.fdfcoefs[n])-2+chan.nonzeroindices[n]), in, chan.fdfcoefs[n], chan.nonzeroindices[n])
        @. rayleighwave += 10^(chan.pathgains[n]/10) * chan.pathrayleighamp[n] * exp(j*2π*(chan.pathdopplers[n]*t-chan.pathoffset[n])) * delayedwave
    end
    return rayleighwave
end

end