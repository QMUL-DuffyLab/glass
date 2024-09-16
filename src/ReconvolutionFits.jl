module ReconvolutionFits

using Interpolations, Symbolics, FFTW, LsqFit, DelimitedFiles, Plots
using Plots.PlotMeasures

"""
given a set of x values and an array p with 2n elements, calculate
an n-exponential decay
"""
function multi_exp(x, p)
    # expects [A1, A2,...τ1, τ2...]
    n_exp = length(p)÷2
    sum([p[i] * exp.(-x .* (1.0 / p[i + n_exp])) for i=1:n_exp])
end

function convol(x, h)
    real.(ifft!(fft(x) .* fft(h)))
end

"""
perform reconvolution fit. X is a 3-column array with time values,
the normalised (unshifted) IRF and the lifetimes (a little more
explanation about this is below - ugly but works). the IRF has to be
interpolated to the time values of the bins first, then shifted, then
convolved with the histogram data
"""
function reconv(X, args)
    # args should be [amplitudes..., irf_shift]
    x = X[:, 1]
    irf = X[:, 2]
    taus = X[:, 3]
    y = zeros(Float64, length(x))
    itp = linear_interpolation(x, irf, extrapolation_bc=Line())
    irf_interp = itp(x .- args[end])
    norm_irf = irf_interp ./ sum(irf_interp)
    for i = 1:length(args) - 1
        y .+= args[i] * exp.(-x ./ taus[i])
    end
    z = convol(y, norm_irf)
end

"""
one function which uses the fit function and fit data to plot
both the simulated trace and corresponding fit for both the tail
and reconvolution fits. the Julia plots backend isn't as fully
featured as matplotlib but it does the job, I think.
"""
function plot_fit(fn, fit, xfit, xdata, ydata, outfile; irf=nothing)
    fontsize = 18
    gr(size=(1000,600), fontfamily="sans-serif", linewidth=3,
       framestyle=:box, label=false, grid=true, tickfontsize=fontsize,
       legend_font_pointsize=fontsize, guidefontsize=fontsize, margin = 5.0mm)
    yfit = fn(xfit, fit.param)
    plot(yscale = :log10, minorgrid = true, xlabel = "time (ns)",
          ylabel = "counts", ylims = (1e-4, 2.0))
    if !isnothing(irf)
        plot!(xdata .- fit.param[end], irf,
              label="IRF", linestyle=:dash, lw=1.5, lc=:black)
    end
    plot!(xdata, ydata, label="data")
    plot!(xdata, yfit, label="fit")
    savefig(outfile)
end

"""
wrapper to do a tail fit and then reconvolution fit for a given
simulated trace. you have to give it a set of starting guesses for the
lifetimes (in nanoseconds); the function will generate the right
number of parameters etc. depending how many lifetimes you give it.
so e.g.
fit(output_file, [1.1])
will do a one-exponential fit, whereas
fit(output_file, [1.1, 0.02])
will do a two-exponential and so on. the IRF file is assumed to be
in the same folder as the output file with the name 'pulse.txt'; if it
isn't, add the IRF filename at the end. but for my simulated traces this
will always be true so don't worry.
Outputs a pair of text files with tail fit and reconvolution fit data,
and corresponding plots showing the simulated data with the fits.
"""
function fit(filename, τᵢ, irf_file=nothing)
    """
    first load in an output histogram and get the sum of the
    emissive counts; these correspond to the output of the experimental
    setup so they're what we want to fit. Note that what I do is pack
    the bin values, the unnormalised counts and the normalised counts,
    use the unnormalised counts to determine errors for each bin,
    then use the normalised counts for the fit.
    We cut the data off at the maximum count first and fit the tail
    to get time constants for the exponentials, then we hold those
    constant and fit the amplitudes with the IRF below.

    NB: this won't work for an experimental histogram because the
    format of the output won't be the same. I'll write a
    sanitise_counts function to check whether it's a simulated or
    experimental trace if necessary, but this works fine if you only
    want to fit my simulated results.
    """
    bcle = readdlm(filename)
    labels = bcle[1, 1:end]
    emissive = bcle[2, 1:end]
    bincounts = bcle[3:end, 1:end]
    # pull the emissive columns and sum row-by-row (bin-by-bin)
    # if there are more than one - otherwise this sum does nothing
    ec = sum(bincounts[:, emissive .> 0], dims=2)
    max_time = bincounts[argmax(ec), 1]
    xyn = hcat(bincounts[:, 1], ec, ec / maximum(ec))
    tail = xyn[xyn[:, 1] .>= max_time, :]
    tail_σ = zeros(Float64, length(tail[:, 2])) 
    for i = 1:length(tail_σ)
        count = tail[i, 2] > 0 ? tail[i, 2] : one(eltype(tail_σ))
        # error for each bin based on count in bin
        tail_σ[i] = sqrt((1.0 / count) + (1.0 / maximum(ec)))
    end
    # then higher weights for lower errors
    wt = 1.0 ./ tail_σ
    n_exp = length(τᵢ)
    # start with evenly weighted amplitudes
    p0 = [[1.0/n_exp for i=1:n_exp]..., τᵢ .* 1e9...]
    # all fitted parameters should be positive
    lbs = [0. for i=1:n_exp*2]
    ubs = [Inf for i=1:n_exp*2]
    x = 1e9 .* (tail[:, 1])
    y = tail[:, 3]

    # do tail fit, get fit data, print out, plot
    tail_fit = curve_fit(multi_exp, x, y, wt, p0, lower=lbs, upper=ubs)
    # try
        cov = estimate_covar(tail_fit)
    # catch
    #     @warn "could not calculate covariances for tail fit"
    #     cov = Inf
    # end
    # try
        fit_error = stderror(tail_fit)
    # catch
    #     @warn "could not calculate error for tail fit"
    #     fit_error = Inf
    # end
    # try
        margin_of_error = margin_error(tail_fit, 0.05)
    # catch
    #     @warn "could not calculate margin of error for tail fit"
    #     margin_of_error = Inf
    # end
    println("tail fit (amplitudes then lifetimes (ns)): $(tail_fit.param)")
    println("covar: $(cov)")
    println("errors: $(fit_error)")
    println("margin: $(margin_of_error)")
    outfile = splitext(filename)[1] * "_tail_fit_n_$(length(τᵢ)).txt"
    open(outfile, "w") do io
        println(io, "best fit (amplitudes then lifetimes (ns)): $(tail_fit.param)")
        println(io, "covar: $(cov)")
        println(io, "errors: $(fit_error)")
        println(io, "errors: $(fit_error)")
        println(io, "margin: $(margin_of_error)")
    end
    outfile = splitext(filename)[1] * "_tail_fit_n_$(length(τᵢ)).pdf"
    println(outfile)
    plot_fit(multi_exp, tail_fit, x, x, y, outfile)

    """
    now time for the reconvolution fit. first, read in the IRF file;
    if no filename's given just pull the pulse output file which will
    be in the same directory as the histogram.
    then we have to wrap up the x values (time bins), normalised IRF
    and our lifetimes into the first argument of the fit function,
    because otherwise they'll be varied by the least-squares fitter.
    the way i do this is to make a long array of zeros and insert the
    fitted lifetimes into that array; this was necessary for the old
    scipy curve_fit script i was using, but i'm not actually sure it is
    now. but it works, either way, although it is extremely ugly.
    then we call the reconv function defined above, which fits the
    amplitudes and IRF shift while holding the lifetimes constant.

    NB: if using an experimental IRF file it should just be a
    two-column text file containing time values and IRF values,
    otherwise readdlm might not work here.
    """
    if !isnothing(irf_file)
        irf_data = readdlm(irf_file)
    else
        irf_file = joinpath(dirname(filename), "pulse.txt")
        irf_data = readdlm(irf_file)
    end
    irf = irf_data[:, 2]
    irf_norm = irf ./ sum(irf)
    itp = linear_interpolation(irf_data[:, 1], irf_norm,
                               extrapolation_bc=Line())
    irf_interp = itp(xyn[:, 1])
    irf_norm = irf_interp ./ maximum(irf_interp)
    taus = zeros(Float64, length(xyn[:, 1]))
    for i = 1:n_exp
        taus[i] = tail_fit.param[i + n_exp]
    end
    irf_shift = 0.0
    max_shift = 100.0
    # wrap the independent variables as described above
    X = hcat(1e9 .* xyn[:, 1], irf_norm, taus)
    # starting parameters are the fitted amplitudes and 0 for IRF shift
    p0 = [tail_fit.param[1:n_exp]..., irf_shift]
    # all fit parameters must be positive
    lbs = [[0.0 for i=1:n_exp]..., -max_shift]
    ubs = [[Inf for i=1:n_exp]..., max_shift]
    reconv_fit = curve_fit(reconv, X, xyn[:, 3], p0, lower=lbs, upper=ubs)
    # try
        cov = estimate_covar(tail_fit)
    # catch
    #     @warn "could not calculate covariances for reconvolution fit"
    #     cov = Inf
    # end
    # try
        fit_error = stderror(tail_fit)
    # catch
    #     @warn "could not calculate error for reconvolution fit"
    #     fit_error = Inf
    # end
    # try
        margin_of_error = margin_error(tail_fit, 0.05)
    # catch
    #     @warn "could not calculate margin of error for reconvolution fit"
    #     margin_of_error = Inf
    # end
    println("reconv fit (amplitudes then IRF shift): $(reconv_fit.param)")
    println("covar: $(cov)")
    println("errors: $(fit_error)")
    println("margin: $(margin_of_error)")
    outfile = splitext(filename)[1] * "_reconv_fit_n_$(length(τᵢ)).txt"
    open(outfile, "w") do io
        println(io, "reconv fit (amplitudes then lifetimes (ns)): $(tail_fit.param)")
        println(io, "covar: $(cov)")
        println(io, "errors: $(fit_error)")
        println(io, "errors: $(fit_error)")
        println(io, "margin: $(margin_of_error)")
    end

    outfile = splitext(filename)[1] * "_reconv_fit_n_$(length(τᵢ)).pdf"
    plot_fit(reconv, reconv_fit, X, xyn[:, 1] .* 1e9, xyn[:, 3],
             outfile, irf=irf_norm)

    # print the fitted trace to a text file
    yreconv = reconv(X, reconv_fit.param)
    outfile = splitext(filename)[1] * "_fit_trace_n_$(length(τᵢ)).txt"
    open(outfile, "w") do io
        writedlm(io, hcat(xyn[:, 1], yreconv))
    end

end

end
