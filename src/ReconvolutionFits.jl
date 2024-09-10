module ReconvolutionFits

using PyPlot, Symbolics, FFTW, LsqFit, DelimitedFiles

function multi_exp(x, p)
    # expects [A1, A2,...τ1, τ2...]
    n_exp = length(p)÷2
    sum([p[i] * exp.(-x .* (1.0 / p[i + n_exp])) for i=1:n_exp])
end

function convol(x, h)
    real.(ifft!(fft(x) .* fft(h)))
end

function reconv(x, irf, τs, amps, irf_shift)

end

function plot_fit(fit, xyn, outfile)
    x = xyn[:, 1]
    x_ns = x .* 1e9
    y = xyn[:, 3]
    yfit = multi_exp(x_ns, fit.param)
    fig, ax = plt.subplots(figsize=(12,8))
    plot(x_ns, y, label="data")
    plot(x_ns, yfit, label="fit")
    ax.set_ylim(1e-5, 2.0) # normalised count max to 1
    ax.set_yscale("log")
    plt.grid(visible=true)
    ax.legend()
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("counts")
    savefig(outfile)
    plt.close()
end

function fit(filename, τᵢ, irf_file, pulse_fwhm, pulse_peak)
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
    fit_res = curve_fit(multi_exp, x, y, wt, p0, lower=lbs, upper=ubs)
    fit_error = stderror(fit_res)
    margin_of_error = margin_error(fit_res, 0.05)

    println("best fit: $(fit_res.param)")
    println("errors: $(fit_error)")
    println("margin: $(margin_of_error)")
    outfile = joinpath(dirname(filename), "fit.pdf")
    println(outfile)
    plot_fit(fit_res, tail, outfile)

end

end
