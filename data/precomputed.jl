"""
    const DIST_1718

A dictionary mapping `(year, month)` pairs to [`SargassumDistribution`](@ref)s for the years 2017 and 2018 \
computed using `SargassumFromAFAI.jl` at months where sufficient data are available.
"""
const DIST_1718 = SargassumDistribution(joinpath(@__DIR__, "dists-2017-2018.nc"))