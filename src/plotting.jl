using Makie, CairoMakie, GeoMakie
using GeoMakie.GeoJSON
using Latexify
########################################################

function default_fig()
    return Figure(
    resolution = (1920, 1080), 
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
            tickformat = values -> [latexify(value, fmt = FancyNumberFormatter(3)) for value in values],
            width = barwidth,
            valign = :top,
            tellheight = false
    )

    rowsize!(data_legend, 2, barlength) # relative height of bar and its label
    # colsize!(fig_pos, 2, Relative(1/10)) # relative width of data legend    
end