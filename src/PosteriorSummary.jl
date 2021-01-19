module PosteriorSummary

export effective_sample_size,
       hpd_interval

import MCMCDiagnostics.effective_sample_size

function hpd_interval(x::AbstractVector{Float64}; alpha::Float64 = 0.05, sorted::Bool = false)
    n = length(x)
    span = Int(round((1.0 - alpha) * n))
    smallest_rng = Inf
    rng_upper = Inf
    rng_lower = -Inf
    sx::AbstractVector{Float64} = sorted ? x : sort(x)
    lower = 0
    for upper = span:n
        lower = upper - span + 1
        rng = sx[upper] - sx[lower]
        if rng < smallest_rng
            smallest_rng = rng
            rng_lower = sx[lower]
            rng_upper = sx[upper]
        end
    end

    return (lower = rng_lower, upper = rng_upper)
end

end
