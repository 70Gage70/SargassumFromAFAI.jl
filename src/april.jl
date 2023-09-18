include(joinpath(@__DIR__, "main.jl"))

file = joinpath(@__DIR__, "..", "data", "afai-2018-04.nc")
ncinfo(file)

params = AFAIParameters()
afai = AFAI(file, params)

coast_mask = CoastMask(afai)
coast_masked!(afai, coast_mask)

afai_median_background = afai_median_img(afai)

classification = pixel_classification(afai, afai_median = afai_median_background);

unmixed = pixel_unmixing(afai, pixel_classification = classification);

lon_bins, lat_bins, coverage_tot = coverage(afai, unmixed = unmixed)
distribution = SargassumDistribution(lon_bins, lat_bins, DateTime(2018, 4), coverage_tot)