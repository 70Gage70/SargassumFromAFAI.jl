module SargassumFromAFAI

# core functionality
using Dates, Statistics, ArgCheck

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
export clean_pacific!, coast_and_clouds!, pixel_classify!, pixel_unmix!
export distribution_to_nc, afai_to_distribution

include(joinpath(@__DIR__, "..", "data", "precomputed.jl"))
export DIST_1718

include(joinpath(@__DIR__, "plotting.jl"))
export coast!, clouds!, plot, plot!

export show # various Base extensions

end # module
