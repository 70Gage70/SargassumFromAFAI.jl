########################################################
# COAST

"""
    coast!(ax, dist; args...)

Add a heatmap of `dist.coast` to `ax::Makie.Axis`, where dist can be [`AFAI`](@ref) or [`SargassumDistribution`](@ref).

### Optional Arguments 

- `args...`: All keyword arguments are passed directly to `Makie.heatmap!`.
"""
function coast!(
    ax::Axis, 
    dist::Union{AFAI, SargassumDistribution}; 
    args...)

    lon = dist.lon
    lat = dist.lat
    coast = dist.coast

    defaults = (
        interpolate = false, 
        colorrange = (0.3, 0.4), 
        lowclip = RGBAf(0,0,0,0),
        highclip = colorant"red"
    )
    
   return heatmap!(ax, lon, lat, coast; merge(defaults, args)...)
end

########################################################
# CLOUDS

"""
    clouds!(ax, dist, week; args...)

Add a heatmap of `dist.clouds[:, :, week]` to `ax::Makie.Axis`, where dist can be [`AFAI`](@ref) or [`SargassumDistribution`](@ref).

### Optional Arguments 

- `args...`: All keyword arguments are passed directly to `Makie.heatmap!`.
"""
function clouds!(
    ax::Axis, 
    dist::Union{AFAI, SargassumDistribution},
    week::Integer; 
    args...)

    @argcheck week in [1, 2, 3, 4]

    lon = dist.lon
    lat = dist.lat
    clouds = dist.clouds[:,:,week]

    defaults = (
        interpolate = false, 
        colorrange = (0.3, 0.4), 
        lowclip = RGBAf(0,0,0,0),
        highclip = colorant"black"
    )
    
   return heatmap!(ax, lon, lat, clouds; merge(defaults, args)...)
end

########################################################
# AFAI

"""
    plot(afai; show_coast, show_clouds)

Plot `afai.afai` for each of the four weeks on one graph.

### Optional Arguments

- `show_coast`: Highlight the coastlines in each graph via [`coast!`](@ref).
- `show_clouds`: Highlight clouds/missing data in each graph via [`clouds!`](@ref).
- `limits`: A `NTuple{4, Int64}` giving the limits of the graph in the form `(lon_min, lon_max, lat_min, lat_max)`. Default `(-100, -38, 0, 35)`.
"""
function plot(
    afai::AFAI; 
    show_coast::Bool = false, 
    show_clouds::Bool = false,
    limits::NTuple{4, Int64} = (-100, -38, 0, 35))

    lon = afai.lon
    lat = afai.lat
    afai_data = afai.afai

    set_theme!(GEO_THEME())
    fig = Figure()

    # Day 7
    ax = Axis(fig[1, 1], title = "Days 1-8", limits = limits)
    heatmap!(ax, lon, lat, afai_data[:,:,1])
    show_coast ? coast!(ax, afai) : nothing
    show_clouds ? clouds!(ax, afai, 1) : nothing
    land!(ax)
    
    # Day 14
    ax = Axis(fig[1, 2], title = "Days 8-15", limits = limits)
    heatmap!(ax, lon, lat, afai_data[:,:,2])
    show_coast ? coast!(ax, afai) : nothing
    show_clouds ? clouds!(ax, afai, 2) : nothing
    land!(ax)
    
    # Day 21
    ax = Axis(fig[2, 1], title = "Days 15-22", limits = limits)
    heatmap!(ax, lon, lat, afai_data[:,:,3])
    show_coast ? coast!(ax, afai) : nothing
    show_clouds ? clouds!(ax, afai, 3) : nothing
    land!(ax)
    
    # Day 28
    ax = Axis(fig[2, 2], title = "Days 22-29", limits = limits)
    heatmap!(ax, lon, lat, afai_data[:,:,4])
    show_coast ? coast!(ax, afai) : nothing
    show_clouds ? clouds!(ax, afai, 4) : nothing
    land!(ax)
    
    return fig
end

########################################################
# SARGASSUM DISTRIBUTION

"""
    plot(sargassum_distribution; show_coast, show_clouds, limits)

Plot `sargassum_distribution` for each of the four weeks in one graph.

### Optional Arguments

- `show_coast`: Highlight the coastlines in each graph via [`coast!`](@ref). Default `false`.
- `show_clouds`: Highlight clouds/missing data in each graph via [`clouds!`](@ref). Default `false`.
- `limits`: A `NTuple{4, Int64}` giving the limits of the graph in the form `(lon_min, lon_max, lat_min, lat_max)`. Default `(-90, -38, -5, 22)`.
"""
function plot(
    sargassum_distribution::SargassumDistribution;
    show_coast::Bool = false,
    show_clouds::Bool = false,
    limits::NTuple{4, Int64} = (-90, -38, -5, 22))

    lon = sargassum_distribution.lon
    lat = sargassum_distribution.lat
    sarg = sargassum_distribution.sargassum

    year_string = Year(sargassum_distribution.time).value
    month_string = monthname(sargassum_distribution.time)

    set_theme!(GEO_THEME())
    fig = Figure(
        size = (1920, 1080),
        fontsize = 50,
        figure_padding = (5, 5, 5, 5))

    # Day 7
    ax = Axis(fig[1, 1], title = "Days 1-8", limits = limits)
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg[:,:,1]))
    heatmap!(ax, lon, lat, sarg[:,:,1], 
        colormap = EUREKA,
        colorrange = sarg_limits,
        lowclip = :white)
    show_coast ? coast!(ax, sargassum_distribution) : nothing
    show_clouds ? clouds!(ax, sargassum_distribution, 1) : nothing
    land!(ax)
    
    # Day 14
    ax = Axis(fig[1, 2], title = "Day 8-15", limits = limits)
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg[:,:,2]))
    heatmap!(ax, lon, lat, sarg[:,:,2], 
        colormap = EUREKA,
        colorrange = sarg_limits,
        lowclip = :white)
    show_coast ? coast!(ax, sargassum_distribution) : nothing
    show_clouds ? clouds!(ax, sargassum_distribution, 2) : nothing
    land!(ax)
    
    # Day 21
    ax = Axis(fig[2, 1], title = "Day 15-22", limits = limits)
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg[:,:,3]))
    heatmap!(ax, lon, lat, sarg[:,:,3], 
        colormap = EUREKA,
        colorrange = sarg_limits,
        lowclip = :white)
    show_coast ? coast!(ax, sargassum_distribution) : nothing
    show_clouds ? clouds!(ax, sargassum_distribution, 3) : nothing
    land!(ax)
    
    # Day 28
    ax = Axis(fig[2, 2], title = "Day 22-29", limits = limits)
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg[:,:,4]))
    heatmap!(ax, lon, lat, sarg[:,:,4], 
        colormap = EUREKA,
        colorrange = sarg_limits,
        lowclip = :white)
    show_coast ? coast!(ax, sargassum_distribution) : nothing
    show_clouds ? clouds!(ax, sargassum_distribution, 4) : nothing
    land!(ax)

    fig[0, :] = Label(fig, "$(month_string) $(year_string)")

    return fig
end

"""
    plot(sargassum_distribution, week; show_coast, show_clouds, limits)

Plot `sargassum_distribution` for the week `week`.

### Optional Arguments

- `show_coast`: Highlight the coastlines in each graph via [`coast!`](@ref). Default `false`.
- `show_clouds`: Highlight clouds/missing data in each graph via [`clouds!`](@ref). Default `false`.
- `limits`: A `NTuple{4, Int64}` giving the limits of the graph in the form `(lon_min, lon_max, lat_min, lat_max)`. Default `(-90, -38, -5, 22)`.
"""
function plot(
    sargassum_distribution::SargassumDistribution,
    week::Integer;
    show_coast::Bool = false,
    show_clouds::Bool = false,
    limits::NTuple{4, Int64} = (-90, -38, -5, 22))

    @argcheck week ∈ [1, 2, 3, 4]

    lon = sargassum_distribution.lon
    lat = sargassum_distribution.lat
    sarg = sargassum_distribution.sargassum[:,:,week]

    year_string = Year(sargassum_distribution.time).value
    month_string = monthname(sargassum_distribution.time)

    set_theme!(GEO_THEME())
    fig = Figure(
        size = (1920, 800),
        fontsize = 50,
        figure_padding = (5, 100, 5, 5))

    ax = Axis(fig[1, 1], title = "$(month_string) $(year_string), week  $(week)", limits = limits)
    
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg))

    heatmap!(ax, lon, lat, sarg, 
        colormap = EUREKA,
        colorrange = sarg_limits,
        lowclip = :white)

    show_coast ? coast!(ax, sargassum_distribution) : nothing
    show_clouds ? clouds!(ax, sargassum_distribution, week) : nothing

    land!(ax)
    
    return fig
end

"""
    sarg!(axis, sargassum_distribution, week; log_scale, args...)

Add a plot of `sargassum_distribution` for the week `week` to `Axis::Makie.Axis`. Returns a `Makie.heatmap!`.

### Optional Arguments

- `log_scale`: Plot on a `log10` scale. Default `false`.
- `args...`: All keyword arguments are passed directly to `Makie.heatmap!`.
"""
function sarg!(
    axis::Axis,
    sargassum_distribution::SargassumDistribution,
    week::Integer;
    log_scale::Bool = false,
    args...)

    @argcheck week ∈ [1, 2, 3, 4]

    lon = sargassum_distribution.lon
    lat = sargassum_distribution.lat
    sarg = sargassum_distribution.sargassum[:,:,week]
    
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg))

    defaults = (
        colormap = EUREKA,
        colorrange = sarg_limits,
        colorscale = log_scale ? log10 : x -> x,
        lowclip = :white
    )

    return heatmap!(axis, lon, lat, sarg; merge(defaults, args)...)
end