using GeoDatases
using Interpolations

lon, lat, land = GeoDatasets.landseamask(resolution = 'c', grid = 5)
land[land .== 2] .= 1 # lake is not ocean, so it's land

itp = GriddedField([:lon, :lat], Dict(:lon => lon, :lat => lat), Dict(:land => land), nothing, nothing, nothing, ref_itp)
itp = itp |> sph2xy |> x -> interpolate(x, interpolant_type = "nearest")