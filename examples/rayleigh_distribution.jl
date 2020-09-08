using MultipathChannel
using Plots
using Random
using StatsBase
using LaTeXStrings
using Statistics

Random.seed!(1)

# Rayleigh Distribution
# ---------- Mathmatical Model ----------
σ = 1/√2
normDistRange = -4σ:0.1:4σ
x = normDistRange
normDist = 1/√(2π*σ^2) * exp.(-x.^2/(2σ^2))

rayleighDistRange = 0:0.05:5σ
x = rayleighDistRange
rayleighDist = x./σ^2 .* exp.(-x.^2/(2σ^2))

# ---------- Rayleigh Fading ----------
channelmodel = "rayleigh"
maxdoppler = 0
samplingrate = 15.36e6
numiters = 50000
amp = Vector{Complex{Float64}}(undef, numiters)
for i = 1:numiters
    chan = initChan(maxdoppler, samplingrate, channelmodel)
    amp[i] = applyChan([1], chan, initialtime=rand()*10000)[1]
end

# Create Histograms of in-phase, quadrature, and envelope
HistoRxI = fit(Histogram, real.(amp), normDistRange)
HistoRxINorm = StatsBase.normalize(HistoRxI, mode=:pdf)
HistoRxQ = fit(Histogram, imag.(amp), normDistRange)
HistoRxQNorm = StatsBase.normalize(HistoRxQ, mode=:pdf)
HistoRxEnv = fit(Histogram, abs.(amp), rayleighDistRange)
HistoRxEnvNorm = StatsBase.normalize(HistoRxEnv, mode=:pdf)

# Show a graph of Histograms
p1 = plot(HistoRxINorm, title="In-phase", legend=false); plot!(normDistRange, normDist); xlabel!("Amplitude"); ylabel!("Power Distribution Function")
p2 = plot(HistoRxQNorm, title="Quadrature", legend=false); plot!(normDistRange, normDist); xlabel!("Amplitude"); ylabel!("Power Distribution Function")
p3 = plot(HistoRxEnvNorm, title="Envelope", legend=false); plot!(rayleighDistRange, rayleighDist); xlabel!("Amplitude"); ylabel!("Power Distribution Function")
plot(p1,p2,p3, layout=(1,3))
savefig("RayleighAmplitude.png")
