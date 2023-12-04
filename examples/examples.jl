"""
    const EXAMPLE_DIST_APRIL_2018

The path to a NetCDF file giving the April 2018 
Sargassum distribution computed using `SargassumFromAFAI.jl`
"""
const EXAMPLE_DIST_APRIL_2018 = joinpath(@__DIR__, "dist-2018-4.nc")

"""
    const EXAMPLE_DIST_MAY_2018

The path to a NetCDF file giving the May 2018 
Sargassum distribution computed using `SargassumFromAFAI.jl`
"""
const EXAMPLE_DIST_MAY_2018 = joinpath(@__DIR__, "dist-2018-5.nc")

"""
    const DIST_2018

The path to a NetCDF file giving the January-July 2018 
Sargassum distributions computed using `SargassumFromAFAI.jl`
"""
const DIST_2018 = SargassumDistribution(joinpath(@__DIR__, "..", "data", "dist-2018.nc"))