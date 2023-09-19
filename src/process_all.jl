using JLD2

include(joinpath(@__DIR__, "main.jl"))

year = 2018
months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]

files = String[]

for i = 1:length(months)
    push!(files, joinpath(@__DIR__, "..", "data", "afai-$(year)-$(months[i]).nc"))
end

dists = []

for file in files
    @info file
    @time dist = afai_to_distribution(file)
    push!(dists, dist)
end

jldsave("dists.jld2", dists = dists)
