########################################################
# GENERAL 

"""
    default_fig()

Create a `Makie.Figure` with a size of `(1920, 1080)`, a fontsize of `50` and a padding 
of `(5, 100, 5, 5)`.
"""
function default_fig()
    return Figure(
    size = (1920, 1080), 
    fontsize = 50,
    figure_padding = (5, 100, 5, 5));
end

"""
    geo_axis(fig_pos; title, limits, xticks, yticks)

Create a `Makie.Axis` suitable for plotting on an equirectangular projection.

### Arguments

- `fig_pos`: A `Makie.GridPosition` where the plot should go. For example if `fig` is a `Makie.Figure`, \
then `fig_pos[1, 1]` puts the axis in the first row and first column of `fig`.

### Optional Arguments 

- `title`: An `AbstractString`. Default `L"\\mathrm{Title}"`, where `Makie.L"..."` creates a `LaTeXString`.
- `limits`: `An NTuple{4, <:Real}` of the form `(xmin, xmax, ymin, ymax)`.
- `xticks`: A list of x tick mark locations.
- `yticks`: A list of y tick mark locations.
"""
function geo_axis(
    fig_pos::GridPosition;
    title::AbstractString = L"\mathrm{Title}",
    limits::Tuple = (-100, -50, 5, 35),
    xticks::Vector{<:Real} = [limits[1], limits[2]],
    yticks::Vector{<:Real} = [limits[3], limits[4]]
    )

    Axis(
        fig_pos,
        limits = limits, 
        title = title,
        xticks = xticks,
        yticks = yticks,
        xticklabelsize = 40,
        yticklabelsize = 40,
        xtickformat = values -> [
            if value > 0 
                L"%$(abs(value)) \, \degree \mathrm{E}" 
            elseif value == 0 
                L"0\degree"
            elseif value < 0
                L"%$(abs(value)) \, \degree \mathrm{W}" 
            end
        for value in values],
        ytickformat = values -> [
            if value > 0 
                L"%$(abs(value)) \, \degree \mathrm{N}" 
            elseif value == 0 
                L"0\degree"
            elseif value < 0
                L"%$(abs(value)) \, \degree \mathrm{S}" 
            end
        for value in values]
    )
end


"""
    land!(axis; landpath, color)

Add a land polygon to `axis::Makie.Axis`. This will be placed on top of any graphics that 
are already on the axis.

### Optional Arguments

- `landpath`: A `String` pointing to the location of the `.geojson` file containing the land features. By \
default this is provided by a Natural Earth 50 m file \
"https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne\\_50m\\_land.geojson" \
which is included automatically. 
- `color`: A `Makie.Colorant` which gives the color of the land. Default is grey via `RGBf(0.5,0.5,0.5)`.
"""
function land!(
    axis::Axis; 
    landpath::String = joinpath(@__DIR__, "..", "geojson", "ne_50m_land.geojson"), 
    color::Colorant = RGBf(0.5,0.5,0.5)) # aka colorant"gray50"

    landpoly = GeoJSON.read(read(landpath, String))
    poly!(axis, landpoly, color = color)
end

"""
    data_legend!(fig_pos, label; colormap, label_fontsize, tick_fontsize, ticks, barlength, barwidth)

Add a `Makie.Colorbar` with label `label` to the `GridPosion` in `fig_pos.`

### Optional Arguments

- `colormap`: The colormap used in the colorbar. Default `SargassumColors.EUREKA`.
- `label_fontsize`: The font size of the label. Default 40.
- `tick_fontsize`: The font size of the label ticks. Default 40.
- `ticks`: A `Vector` of tick values. Default `[0.0, 0.2, 0.4, 0.6, 0.8, 1.0]`.
- `barlength`: The length of the colorbar as a `Relative` size of the grid. Default `Relative(9/10).`
- `barwidth`: The width of the colorbar as a `Relative` size of the grid. Default `Relative(3/10).`
"""
function data_legend!(
    fig_pos::GridPosition,
    label::AbstractString = L"% \, \mathrm{Covr.}";
    colormap = EUREKA,
    label_fontsize::Real = 40,
    tick_fontsize::Real = 40,
    ticks::Vector{<:Real} = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
    barlength::Relative = Relative(9/10),
    barwidth::Relative = Relative(3/20))

    data_legend = GridLayout(fig_pos)

    Label(
        data_legend[1,1], 
        label,
        # L"\left(\frac{\pi_{\text{stat}}}{%$(latexify(colors_max, fmt = FancyNumberFormatter(3), env = :raw))}\right)^{1/4}", 
        fontsize = label_fontsize,
        valign = :bottom, 
        tellheight = false
    )
    
    Colorbar(
            data_legend[2, 1], 
            limits = extrema(ticks),
            colormap = colormap,
            ticklabelsize = tick_fontsize, 
            ticks = ticks, 
            tickformat = values -> [latexify(100*value, fmt = FancyNumberFormatter(3)) for value in values],
            width = barwidth,
            valign = :top,
            tellheight = false
    )

    rowsize!(data_legend, 2, barlength) # relative height of bar and its label
    # colsize!(fig_pos, 2, Relative(1/10)) # relative width of data legend    
end

########################################################
# COAST MASK

"""
    plot(coast_mask)

Plot the coastlines in `coast_mask`.
"""
function plot(coast_mask::CoastMask)
    lon = coast_mask.lon
    lat = coast_mask.lat
    mask = coast_mask.mask

    fig = default_fig()

    ax = geo_axis(fig[1, 1], title = L"\text{Coastlines}", limits = (-100, -38, 0, 35))
    
    heatmap!(ax, lon, lat, mask, 
        interpolate = false, 
        nan_color = :black, 
        colorrange = (0.3, 0.4), 
        highclip = :gray)
    land!(ax)
    
   return fig
end

########################################################
# AFAI

"""
    plot(afai; coast_mask)

Plot `afai.afai` for each of the four weeks on one graph.

### Optional Arguments

- `coast_mask`: If provided, the [`coast_masked`](@ref) version of `afai.afai` is plotted. Default `nothing`.
"""
function plot(afai::AFAI; coast_mask::Union{Nothing, CoastMask} = nothing)
    lon = afai.lon
    lat = afai.lat

    if coast_mask === nothing
        afai_data = afai.afai
    else
        afai_data = coast_masked(afai, coast_mask)
    end

    fig = default_fig()

    # Day 8
    ax = geo_axis(fig[1, 1], title = L"\text{Days 1-8}", limits = (-100, -38, 0, 35))
    heatmap!(ax, lon, lat, afai_data[:,:,1])
    land!(ax)
    
    # Day 15
    ax = geo_axis(fig[1, 2], title = L"\text{Day 9-15}", limits = (-100, -38, 0, 35))
    heatmap!(ax, lon, lat, afai_data[:,:,2])
    land!(ax)
    
    # Day 22
    ax = geo_axis(fig[2, 1], title = L"\text{Day 16-22}", limits = (-100, -38, 0, 35))
    heatmap!(ax, lon, lat, afai_data[:,:,3])
    land!(ax)
    
    # Day 29
    ax = geo_axis(fig[2, 2], title = L"\text{Day 23-29}", limits = (-100, -38, 0, 35))
    heatmap!(ax, lon, lat, afai_data[:,:,4])
    land!(ax)
    
    fig
end

########################################################
# SARGASSUM DISTRIBUTION

"""
    plot(sargassum_distribution; limits, size, legend)

Plot `sargassum_distribution` for each of the four weeks in one graph.

### Optional Arguments

- `limits`: A `NTuple{4, Int64}` giving the limits of the graph in the form `(lon_min, lon_max, lat_min, lat_max)`. Default `(-90, -38, -5, 22)`.
- `size`: A `NTuple{2, Int64}` giving the size of the figure. Default `(1920, 1080)`.
- `legend`: A `Bool`, displays the legends if `true`. Default `true`.
"""
function plot(
    sargassum_distribution::SargassumDistribution;
    limits::NTuple{4, Int64} = (-90, -38, -5, 22),
    size::NTuple{2, Int64} = (1920, 1080),
    legend::Bool = true)

    lon = sargassum_distribution.lon
    lat = sargassum_distribution.lat
    sarg = sargassum_distribution.sargassum

    year_string = Year(sargassum_distribution.time).value
    month_string = monthname(sargassum_distribution.time)

    fig = Figure(
        size = size,
        fontsize = 50,
        figure_padding = (5, 5, 5, 5))

    # Day 8
    ax = geo_axis(fig[1, 1], title = L"\text{Days 1-8}", limits = limits)
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg[:,:,1]))
    heatmap!(ax, lon, lat, sarg[:,:,1], 
        colormap = EUREKA,
        colorrange = sarg_limits,
        lowclip = :white)
    land!(ax)
    
    # Day 15
    ax = geo_axis(fig[1, 2], title = L"\text{Day 9-15}", limits = limits)
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg[:,:,2]))
    heatmap!(ax, lon, lat, sarg[:,:,2], 
        colormap = EUREKA,
        colorrange = sarg_limits,
        lowclip = :white)
    land!(ax)
    
    # Day 22
    ax = geo_axis(fig[2, 1], title = L"\text{Day 16-22}", limits = limits)
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg[:,:,3]))
    heatmap!(ax, lon, lat, sarg[:,:,3], 
        colormap = EUREKA,
        colorrange = sarg_limits,
        lowclip = :white)
    land!(ax)
    
    # Day 29
    ax = geo_axis(fig[2, 2], title = L"\text{Day 23-29}", limits = limits)
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg[:,:,4]))
    heatmap!(ax, lon, lat, sarg[:,:,4], 
        colormap = EUREKA,
        colorrange = sarg_limits,
        lowclip = :white)
    land!(ax)

    if legend
        sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg))
        data_legend!(fig[:, 3], 
            ticks = collect(range(sarg_limits[1], sarg_limits[2], length = 4))
        )

        colsize!(fig.layout, 3, Relative(2/15)) # relative width of data legend 
    end

    fig[0, :] = Label(fig, L"\text{%$(month_string) %$(year_string)}")

    return fig
end

"""
    plot(sargassum_distribution, week; limits, size, legend)

Plot `sargassum_distribution` for the week `week`.

### Optional Arguments

- `limits`: A `NTuple{4, Int64}` giving the limits of the graph in the form `(lon_min, lon_max, lat_min, lat_max)`. Default `(-90, -38, -5, 22)`.
- `size`: A `NTuple{2, Int64}` giving the size of the figure. Default `(1920, 1080)`.
- `legend`: A `Bool`, displays the legends if `true`. Default `true`.
"""
function plot(
    sargassum_distribution::SargassumDistribution,
    week::Integer;
    title::Union{Nothing, AbstractString} = nothing,
    limits::NTuple{4, Int64} = (-90, -38, -5, 22),
    size::NTuple{2, Int64} = (1920, 800),
    legend::Bool = true)

    @assert week ∈ [1, 2, 3, 4]

    lon = sargassum_distribution.lon
    lat = sargassum_distribution.lat
    sarg = sargassum_distribution.sargassum[:,:,week]

    year_string = Year(sargassum_distribution.time).value
    month_string = monthname(sargassum_distribution.time)

    fig = Figure(
        size = size,
        fontsize = 50,
        figure_padding = (5, 100, 5, 5))

    if title === nothing
        tit = L"\text{%$(month_string) %$(year_string)}, \, \mathrm{week} \, %$(week)"
    else 
        tit = title
    end

    ax = geo_axis(fig[1, 1], title = tit, limits = limits)
    
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg))

    heatmap!(ax, lon, lat, sarg, 
        colormap = EUREKA,
        colorrange = sarg_limits,
        lowclip = :white)

    if legend
        data_legend!(fig[1, 2], 
            ticks = collect(range(sarg_limits[1], sarg_limits[2], length = 4))
        )

        colsize!(fig.layout, 2, Relative(1/15)) # relative width of data legend 
    end

    land!(ax)
    
    return fig
end

"""
    plot!(axis. sargassum_distribution, week)

Add a plot of `sargassum_distribution` for the week `week` to `Axis::Makie.Axis`.
"""
function plot!(
    axis::Axis,
    sargassum_distribution::SargassumDistribution,
    week::Integer)

    @assert week ∈ [1, 2, 3, 4]

    lon = sargassum_distribution.lon
    lat = sargassum_distribution.lat
    sarg = sargassum_distribution.sargassum[:,:,week]
    
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg))

    return heatmap!(axis, lon, lat, sarg, 
        colormap = EUREKA,
        colorrange = sarg_limits,
        colorscale = log10,
        lowclip = :white)
end