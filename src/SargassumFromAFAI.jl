module SargassumFromAFAI

export AFAIParameters, AFAI, CoastMask, SargassumDistribution
export coast_masked, coast_masked!, afai_median, pixel_classification, pixel_unmixing, coverage, afai_to_distribution
export plot
export download_data, data_path, data_rm

include(joinpath(@__DIR__, "main.jl"))
include(joinpath(@__DIR__, "..", "data", "remote-files.jl"))

end # module
