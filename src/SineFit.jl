module SineFit

import LsqFit
import FFTW
import Statistics

export estimate_wave_parameters, calculate_wave_shape, sin

struct WaveFitParams
    amplitude::Float64
    frequency::Float64
    phaseshft::Float64
    vertoffst::Float64
end

"""
    estimate_wave_parameters(xdata::Vector{T}, ydate::Vector{T};
        sample_rate=nothing)::WaveFitParams where {T<:Number}

Use Fourier analysis and statistics to find the best monochromatic sine wave parameters for data.
`xsym` and `ysym` are the column symbols for the x and y columns within DataFrame `df`, respectively.
Returns a `WaveFitParams` struct with frequency and phase shift as circular frequencies resp. radians.

If `sample_rate` is not specified, the rate is determined automatically from the `xdata` argument.
"""
function estimate_wave_parameters(xdata::Vector{T}, ydata::Vector{T};
        sample_rate=nothing)::WaveFitParams where {T<:Number}
    if sample_rate === nothing
        sample_rate = 1 / Statistics.mean(xdata[2:end] - xdata[1:end-1])
    end

    offset = (maximum(ydata)+minimum(ydata))/2
    ampl = Statistics.mean(Statistics.quantile(abs.(ydata .- offset), 99/100))
    fft = FFTW.rfft(ydata .- offset)
    fftfreq = FFTW.rfftfreq(length(ydata), sample_rate)

    # Determine dominant frequency and calculate phase shift.
    maxfreqloc = argmax(abs.(fft))
    maxfreq = fftfreq[maxfreqloc]
    shift = mod(atan(imag(fft[maxfreqloc]), real(fft[maxfreqloc])) + pi/2, 2pi)
    freq = maxfreq * 2pi

    WaveFitParams(ampl, freq, shift, offset)
end

function sin_model(xs, params)
    params[1] .* Base.sin.(params[2] .* xs .+ params[3]) .+ params[4]
end

"""
    calculate_wave_shape(xdata::Vector{T}, ydate::Vector{T};
        sample_rate=nothing)::WaveFitParams where {T<:Number}

Fit a sine function (`sin_model`, i.e. `a sin(b x + c) + d`) to the given data `xdata`/`ydata`.

This function first estimates the parameters using normal statistics and fourier analysis using `estimate_wave_parameters()`.
Then, a least-squares fit (`LsqFit`) is used to refine the parameters -- for example, to adjust the frequency, as a FFT has
a very limited frequency resolution. This helps a lot with accuracy.

If `sample_rate` is not specified, the rate is determined automatically from the `xdata` argument.
"""
function calculate_wave_shape(xdata::Vector{T}, ydata::Vector{T};
        sample_rate=nothing)::WaveFitParams where {T<:Number}

    estimate = estimate_wave_parameters(xdata, ydata, sample_rate=sample_rate)
    initial = [estimate.amplitude, estimate.frequency, estimate.phaseshft, estimate.vertoffst]

    fit = LsqFit.curve_fit(sin_model, xdata, ydata, initial)
    params = fit.param

    if params[1] < 0
        params[1] *= -1
        params[3] -= 180
    end
    params[3] = mod(params[3], 360)
    if params[3] > 180
        params[3] = params[3] - 360
    elseif params[3] < -180
        params[3] = params[3] + 360
    end

    WaveFitParams(params...)
end


"""
    sin(x::Vector{T}, wp::WaveFitParams)::Vector{T} where {T<:Number}

Return sin for `x` with the parameters given in `wp`. Useful in order to plot the fit function, for example.

```
wp = calculate_wave_shape(xs, ys)
Plots.plot(xs, ys)
Plots.plot!(xs, sin.(xs, wp))
```

"""
function sin(xs::Vector{T}, wp::WaveFitParams)::Vector{T} where {T<:Number}
    wp.amplitude .* Base.sin.(wp.frequency .* xs .+ wp.phaseshft) .+ wp.vertoffst
end

end # module
