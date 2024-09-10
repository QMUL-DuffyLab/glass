module ReconvolutionFits

using PythonPlot, Interpolations, Symbolics, FFTW, LsqFit, DelimitedFiles

function multi_exp(x, p)
    # expects [A1, A2,...τ1, τ2...]
    n_exp = length(p)÷2
    sum([p[i] * exp.(-x .* (1.0 / p[i + n_exp])) for i=1:n_exp])
end

function convol(x, h)
    real.(ifft!(fft(x) .* fft(h)))
end

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

function plot_fit(fn, fit, xfit, xdata, ydata, outfile; irf=nothing)
    x_ns = xdata .* 1e9
    yfit = fn(xfit, fit.param)
    fig, ax = pyplot.subplots(figsize=(12,8))
    plot(x_ns, ydata, label="data")
    plot(x_ns, yfit, label="fit")
    if !isnothing(irf)
        plot(x_ns .- fit.param[end], irf, label="IRF")
    end
    ax.set_ylim(1e-5, 2.0) # normalised count max to 1
    ax.set_yscale("log")
    pyplot.grid(visible=true)
    ax.legend()
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("counts")
    savefig(outfile)
    plotclose()
end

function fit(filename, τᵢ, irf_file)
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
        tail_σ[i] = sqrt((1.0 / count) + (1.0 / maximum(ec)))
    end
    wt = 1.0 ./ tail_σ
    n_exp = length(τᵢ)
    p0 = [[1.0/n_exp for i=1:n_exp]..., τᵢ .* 1e9...]
    lbs = [0. for i=1:n_exp*2]
    ubs = [Inf for i=1:n_exp*2]
    # x = 1e9 .* (tail[:, 1] .- max_time)
    x = 1e9 .* (tail[:, 1])
    y = tail[:, 3]
    tail_fit = curve_fit(multi_exp, x, y, wt, p0, lower=lbs, upper=ubs)
    cov = estimate_covar(tail_fit)
    fit_error = stderror(tail_fit)
    margin_of_error = margin_error(tail_fit, 0.05)
    println("tail fit (amplitudes then lifetimes (ns)): $(tail_fit.param)")

    outfile = joinpath(dirname(filename), "tail_fit.txt")
    open(outfile, "w") do io
        println(io, "best fit (amplitudes then lifetimes (ns)): $(tail_fit.param)")
        println(io, "covar: $(cov)")
        println(io, "errors: $(fit_error)")
        println(io, "errors: $(fit_error)")
        println(io, "margin: $(margin_of_error)")
    end

    outfile = joinpath(dirname(filename), "tail_fit.pdf")
    println(outfile)
    plot_fit(multi_exp, tail_fit, x, x, y, outfile)

    """
    explain
    """
    irf_data = readdlm(irf_file)
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
    X = hcat(1e9 .* xyn[:, 1], irf_norm, taus)
    p0 = [tail_fit.param[1:n_exp]..., irf_shift]
    lbs = [[0.0 for i=1:n_exp]..., -max_shift]
    ubs = [[Inf for i=1:n_exp]..., max_shift]
    reconv_fit = curve_fit(reconv, X, xyn[:, 3], p0, lower=lbs, upper=ubs)
    # fit_error = stderror(reconv_fit)
    # margin_of_error = margin_error(reconv_fit, 0.05)
    #
    println("reconv fit (amplitudes then IRF shift): $(reconv_fit.param)")
    # println("errors: $(fit_error)")
    # println("margin: $(margin_of_error)")
    outfile = joinpath(dirname(filename), "reconv_fit.txt")
    open(outfile, "w") do io
        println(io, "reconv fit (amplitudes then lifetimes (ns)): $(tail_fit.param)")
        # println(io, "covar: $(cov)")
        # println(io, "errors: $(fit_error)")
        # println(io, "errors: $(fit_error)")
        # println(io, "margin: $(margin_of_error)")
    end

    outfile = joinpath(dirname(filename), "reconv_fit.pdf")
    plot_fit(reconv, reconv_fit, X, xyn[:, 1], xyn[:, 3],
             outfile, irf=irf_norm)
end

end
