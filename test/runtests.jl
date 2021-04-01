using Test
using SineFit

function wp_roughly_equal(a, b)
    ok = true
    for f in [:amplitude, :frequency, :phaseshft, :vertoffst]
        ok &= round(getfield(a, f), digits=3) == round(getfield(b, f), digits=3)
    end
    ok
end

function test_basic_estimation()
    step = 0.01
    xs = collect(range(0, 10, step=step))
    ys = 3 .* Base.sin.(2pi * 1 .* xs .+ 2.2) .- 5
    wp = SineFit.estimate_wave_parameters(xs, ys, sample_rate=round(Int32, 1/step))

    println("test_basic_estimation: ", wp)
    wp
end

@test wp_roughly_equal(test_basic_estimation(), SineFit.WaveFitParams(3.0, 6.276908398780805, 2.230908739161565, -5.0))

function test_basic_fit()
    step = 1.
    xs = collect(range(0, 1000, step=step))
    ys = .3 .* Base.sin.(2pi * .1 .* xs .+ 1.1) .- 5
    wp = SineFit.calculate_wave_shape(xs, ys, sample_rate=round(Int32, 1/step))

    println("test_basic_fit: ", wp)
    wp
end

@test wp_roughly_equal(test_basic_fit(), SineFit.WaveFitParams(0.299999999999996, 0.6283185307179554, 1.1, -5.0))

