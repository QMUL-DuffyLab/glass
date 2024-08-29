module ReconvolutionFits

using Symbolics, FFTW, LsqFit

function multi_exp(t, p)
    # expects [A1, τ1, A2, τ2...]
    sum([p[2*i - 1] * exp.(p[2*i] * t) for i=1:(length(p)÷2)])
end



end
