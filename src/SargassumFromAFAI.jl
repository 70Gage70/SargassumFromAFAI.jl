module SargassumFromAFAI

# core functionality
using Dates, Statistics

# i/o
using NetCDF

# used for creating the coast mask (GeoDatasets.landseamask + interpolation)
using GeoDatasets
using Interpolations

# median filter
using ImageFiltering

# removing pacific datapoints
using PolygonInbounds

# plotting utilities
using SargassumColors
using Makie, CairoMakie, GeoMakie
using GeoMakie.GeoJSON
using Latexify

# downloading and maintaining AFAI data
using RemoteFiles

###############################################################

include(joinpath(@__DIR__, "..", "data", "remote-files.jl"))
export download_data, data_path, data_rm

include(joinpath(@__DIR__, "time.jl"))
export TREF, time2months, months2time

include(joinpath(@__DIR__, "..", "data", "earth-polygons.jl"))
export VERTICES_PACIFIC_PANAMA, VERTICES_NORTH_ATLANTIC

include(joinpath(@__DIR__, "main.jl"))
export AFAIParameters, AFAI, SargassumDistribution
export clean_pacific!, coast_and_clouds!
export coast_masked, coast_masked!, afai_median, pixel_classification, pixel_unmixing, coverage, monthly_total
export distribution_to_nc, afai_to_distribution

include(joinpath(@__DIR__, "..", "examples", "examples.jl"))
export EXAMPLE_DIST_APRIL_2018, EXAMPLE_DIST_MAY_2018, DIST_2018

include(joinpath(@__DIR__, "plotting.jl"))
export plot, plot!

export show # various Base extensions

end # module
