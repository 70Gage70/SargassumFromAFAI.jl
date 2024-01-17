using Makie, CairoMakie, GeoMakie
using GeoMakie.GeoJSON
using Latexify

########################################################
# GENERAL 

function default_fig()
    return Figure(
    size = (1920, 1080), 
    fontsize = 50,
    figure_padding = (5, 100, 5, 5));
end

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

function land!(
    axis::Axis; 
    landpath::String = joinpath(@__DIR__, "..", "geojson", "ne_50m_land.geojson"), 
    color::Colorant = RGBf(0.5,0.5,0.5)) # aka colorant"gray50"

    landpoly = GeoJSON.read(read(landpath, String))
    poly!(axis, landpoly, color = color)
end

function data_legend!(
    fig_pos::GridPosition,
    label::AbstractString = L"% \, \mathrm{Covr.}";
    colormap::Union{Symbol, Reverse{Symbol}} = Reverse(:RdYlGn),
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

function plot(coast_mask::CoastMask)
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


function plot(
    sargassum_distribution::SargassumDistribution;
    limits::NTuple{4, Int64} = (-90, -38, -5, 22),
    resolution::NTuple{2, Int64} = (1920, 1080),
    legend::Bool = true)

    lon = sargassum_distribution.lon
    lat = sargassum_distribution.lat
    sarg = sargassum_distribution.sargassum

    year_string = Year(sargassum_distribution.time).value
    month_string = monthname(sargassum_distribution.time)

    fig = Figure(
        size = resolution,
        fontsize = 50,
        figure_padding = (5, 5, 5, 5))

    # Day 8
    ax = geo_axis(fig[1, 1], title = L"\text{Days 1-8}", limits = limits)
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg[:,:,1]))
    heatmap!(ax, lon, lat, sarg[:,:,1], 
        colormap = Reverse(:RdYlGn),
        colorrange = sarg_limits,
        lowclip = :white)
    land!(ax)
    
    # Day 15
    ax = geo_axis(fig[1, 2], title = L"\text{Day 9-15}", limits = limits)
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg[:,:,2]))
    heatmap!(ax, lon, lat, sarg[:,:,2], 
        colormap = Reverse(:RdYlGn),
        colorrange = sarg_limits,
        lowclip = :white)
    land!(ax)
    
    # Day 22
    ax = geo_axis(fig[2, 1], title = L"\text{Day 16-22}", limits = limits)
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg[:,:,3]))
    heatmap!(ax, lon, lat, sarg[:,:,3], 
        colormap = Reverse(:RdYlGn),
        colorrange = sarg_limits,
        lowclip = :white)
    land!(ax)
    
    # Day 29
    ax = geo_axis(fig[2, 2], title = L"\text{Day 23-29}", limits = limits)
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg[:,:,4]))
    heatmap!(ax, lon, lat, sarg[:,:,4], 
        colormap = Reverse(:RdYlGn),
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

function plot(
    sargassum_distribution::SargassumDistribution,
    week::Integer;
    title::Union{Nothing, AbstractString} = nothing,
    limits::NTuple{4, Int64} = (-90, -38, -5, 22),
    resolution::NTuple{2, Int64} = (1920, 800),
    legend::Bool = true)

    @assert week ∈ [1, 2, 3, 4]

    lon = sargassum_distribution.lon
    lat = sargassum_distribution.lat
    sarg = sargassum_distribution.sargassum[:,:,week]

    year_string = Year(sargassum_distribution.time).value
    month_string = monthname(sargassum_distribution.time)

    fig = Figure(
        size = resolution,
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
        colormap = Reverse(:RdYlGn),
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
        colormap = Reverse(:RdYlGn),
        colorrange = sarg_limits,
        lowclip = :white)
end