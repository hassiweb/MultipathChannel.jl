using MultipathChannel
using Plots
using Random

Random.seed!(1)

# TDL-C channel
samplingrate = 15.36e6
rmsdelay = 300e-9
initialtime = 0.0
maxdoppler = 111.1
channelmodel = "TDLC"
rmsdelay = 300e-9
chan = initChan(maxdoppler, samplingrate, channelmodel, rmsdelay=rmsdelay)

# Impulse
impulse = zeros(Float64, 10000)
impulse[1] = 1.0
response = applyChan(impulse, chan)
plot(1:50, impulse[1:50], line=:stem, legend=false)
plot!(1:50, response[1:50], line=:stem, legend=false)
xlabel!("Delay (in sample)"); ylabel!("Envelope")
savefig("TDLC_impulse.png")

# Square wave
sqwave = repeat([1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0]', 20)
response = applyChan(sqwave, chan)
plot(1:50, sqwave[1:50], legend=false)
plot!(1:50, real(response[1:50]), legend=false)
xlabel!("Delay (in sample)"); ylabel!("Envelope")
savefig("TDLC_square_wave.png")

# Sine wave
t = 0:0.05:50
sinewave = sin.(2π*t)
response = applyChan(sinewave, chan)
plot(1:50, real(sinewave[1:50]), legend=false)
plot!(1:50, real(response[1:50]), legend=false)
xlabel!("Delay (in sample)"); ylabel!("Envelope")
savefig("TDLC_sine_wave.png")
