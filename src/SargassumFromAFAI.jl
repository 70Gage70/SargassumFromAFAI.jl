module SargassumFromAFAI

export AFAIParameters, AFAI, CoastMask, SargassumDistribution
export coast_masked, coast_masked!, afai_median, pixel_classification, pixel_unmixing, coverage, monthly_total
export distribution_to_nc, afai_to_distribution
export plot
export download_data, data_path, data_rm
export TREF, time2months, months2time
export EXAMPLE_DIST_APRIL_2018, EXAMPLE_DIST_MAY_2018

include(joinpath(@__DIR__, "main.jl"))
include(joinpath(@__DIR__, "..", "data", "remote-files.jl"))
include(joinpath(@__DIR__, "..", "examples", "examples.jl"))

end # module
