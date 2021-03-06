using MultipathChannel
using Test
using Statistics

@testset "MultipathChannel.jl" begin

    @testset "Average power after fading" begin
        maxdoppler = 0
        samplingrate = 15.36e6
        channelmodel = "TDLC"
        rmsdelay = 300e-9

        t = 0:0.05:50
        sinewave = sin.(2π*t)

        numiters = 1000
        power = Vector{Float64}(undef, numiters)
        for i = 1:numiters
            chan = initChan(maxdoppler, samplingrate, channelmodel, rmsdelay=rmsdelay)
            response = applyChan(sinewave, chan)
            power[i] = sum((abs.(response)).^2)/length(sinewave)
        end
        @test abs(mean(power) - mean(sum((abs.(sinewave)).^2)/length(sinewave))) < 0.05
    end

    @testset "Expected value of Rayleigh distribution" begin
        maxdoppler = 0
        samplingrate = 15.36e6
        channelmodel = "rayleigh"

        numiters = 1000
        amp = Vector{Complex{Float64}}(undef, numiters)
        for i = 1:numiters
            chan = initChan(maxdoppler, samplingrate, channelmodel)
            amp[i] = applyChan([1], chan, initialtime=rand()*10000)[1]
        end
        σ = 1/√2
        expected = σ*√(π/2)
        @test abs(mean(abs.(amp)) - expected) < 0.05
    end

end
