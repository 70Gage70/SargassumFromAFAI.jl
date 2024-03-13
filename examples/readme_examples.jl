# run this code to generate the figures in the readme

using SargassumFromAFAI
using Makie, CairoMakie

dist_april_2018 = DIST_1718[(2018, 4)]
fig = SargassumFromAFAI.plot(dist_april_2018)
rm("figs/april-2018-weeks.png", force = true)
save("figs/april-2018-weeks.png", fig)

fig = SargassumFromAFAI.plot(dist_april_2018, 1)
rm("figs/april-2018-week-1.png", force = true)
save("figs/april-2018-week-1.png", fig)

dist_may_2018 = DIST_1718[(2018, 5)]
fig = SargassumFromAFAI.plot(dist_may_2018)
rm("figs/may-2018-weeks.png", force = true)
save("figs/may-2018-weeks.png", fig)

fig = Figure()
ax = Axis(fig[1, 1])
SargassumFromAFAI.plot!(ax, dist_april_2018, 1)
rm("figs/april-2018-custom.png", force = true)
save("figs/april-2018-custom.png", fig)

# longer calculation
new_params = AFAIParameters(distribution_quant = 0.5)
path_may_2018 = data_path(2018, 5)
new_dist = afai_to_distribution(path_may_2018, new_params)
fig = SargassumFromAFAI.plot(new_dist)
rm("figs/may-2018-weeks-params.png", force = true)
save("figs/may-2018-weeks-params.png", fig)

