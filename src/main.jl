using NetCDF
using Dates
using Statistics
using GeoDatasets
using Interpolations

include("plotting.jl")

#####################################################################################################################
#####################################################################################################################

"""
    struct AFAIParameters{U, T}

A container for the parameters required to process the AFAI data.

### Fields 

- `window_size_coast_mask`: 
- `window_size_median_filter`: 
- `threshold_median`: 
- `afai_U0`:
- `afai_L0`:

### Constructors 

Use `AFAIParameters(; params...)` where each field has a named kwarg.
"""
struct AFAIParameters{U<:Integer, T<:Real}
    window_size_coast_mask::U
    window_size_median_filter::U
    threshold_median::T
    afai_U0::T 
    afai_L0::T


    function AFAIParameters(;
        window_size_coast_mask = 20,
        window_size_median_filter = 50,
        threshold_median = 1.79e-4,
        afai_U0 = 4.41e-2,
        afai_L0 = -8.77e-4)

        return new{eltype(window_size_coast_mask), eltype(threshold_median)}(
            window_size_coast_mask,
            window_size_median_filter,
            threshold_median,
            afai_U0,
            afai_L0)
    end
end

"""
    mutable struct AFAI{U, R, T}

A container for the AFAI data.

### Fields 

- `lon`: A `Vector` of longitudes.
- `lat`: A `Vector` of latitudes.
- `time`: A `Vector` of `DateTime`s.
- `afai`: An array of AFAI values of the form `afai[lon, lat, time]`.
- `params`: A `AFAIParameters`.

### Constructor

Use `AFAI(filename::String, params::AFAIParameters)` where `filename` is a NetCDF of the form `name.nc`.

It is assumed that the file is obtained from the NOAA database:

https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_atlantic_oceanwatch_AFAI_7D.html

### Plotting 

Use `plot(afai::AFAI; coast_mask::Union{Nothing, CoastMask} = nothing)`.
"""
mutable struct AFAI{U<:Integer, T<:Real, R<:Real}
    lon::Vector{T}
    lat::Vector{T}
    time::Vector{DateTime}
    afai::Array{R, 3}
    params::AFAIParameters{U, T}

    function AFAI(filename::String, params::AFAIParameters)
        extension = infile[findlast(==('.'), infile)+1:end]
        @assert extension == "nc" "Require a NetCDF (.nc) file."

        tref = DateTime(1970, 1, 1, 0, 0, 0)
        time = tref + Second.(ncread(filename, "time"))
    
        lon = ncread(filename, "longitude")
        lat = ncread(filename, "latitude")
        afai = ncread(filename, "AFAI")
    
        return new{eltype(params.window_size_coast_mask), eltype(lon), eltype(afai)}(lon, lat, time, afai, params)
    end
end

"""
    struct CoastMask{T}

A container for a multiplicative mask defining the coastlines. 

A `CoastMask` contains one field, `mask`, which should have the same size as `AFAI.afai`. 

The mask is multiplicative in the sense that its entries are `1.0` away from the coast 
and `NaN` on the coast such that multiplying `AFAI.afai[:,:,t]` by `CoastMask.mask` replaces
the entries of `AFAI.afai[:,:,t]` on the coast by `NaN` and leaves other entries unchanged.

### Constructors

Use `CoastMask(afai::AFAI)`.

### Convenience Function 

Use `coast_masked(afai::AFAI, coast_mask::CoastMask)` to create and return the coast masked
version of `AFAI.afai`.

### Plotting 

Use plot(coast_mask::CoastMask).
"""
struct CoastMask{T<:Real}
    mask::Array{T,3}
end

function CoastMask(afai::AFAI)
    lon = afai.lon
    lat = afai.lat
    afai_data = afai.afai
    window_size = afai.params.window_size_coast_mask

    lon_lsm, lat_lsm, lsm = GeoDatasets.landseamask(resolution = 'l', grid = 5)
    lsm[lsm .== 2] .= 1 # lake is not ocean, so it's land
    landseamask_itp = scale(Interpolations.interpolate(lsm, BSpline(Constant())), lon_lsm, lat_lsm)
    landseamask_gridded = [landseamask_itp(lon_i, lat_i) for lon_i in lon, lat_i in lat]

    landseamask_gridded_coast = zeros(eltype(afai_data), size(landseamask_gridded))
    
    for lon_i = 1:size(landseamask_gridded, 1)
        for lat_i = 1:size(landseamask_gridded, 2)
            lons = max(1, lon_i - round(Integer, window_size/2)):min(length(lon), lon_i + round(Integer, window_size/2))
            lats = max(1, lat_i - round(Integer, window_size/2)):min(length(lat), lat_i + round(Integer, window_size/2))
            if 1.0 in landseamask_gridded[lons, lats]
                landseamask_gridded_coast[lon_i, lat_i] = NaN
            else
                landseamask_gridded_coast[lon_i, lat_i] = 1.0
            end
        end
    end

    return CoastMask(landseamask_gridded_coast)
end

"""
    coast_masked(afai::AFAI, coast_mask::CoastMask)

Return an array equal to `afai.afai` masked by `coast_mask.mask`.
"""
function coast_masked(afai::AFAI, coast_mask::CoastMask)
    afai_unmasked = afai.afai
    afai_masked = zeros(eltype(afai_unmasked), size(afai_unmasked))
    mask = coast_mask.mask

    for t_i = 1:length(time)
        afai_masked[:,:,t_i] = afai_unmasked[:,:,t_i] .* mask
    end

    return afai_masked
end

"""
    coast_masked!(afai::AFAI, coast_mask::CoastMask)

Update `afai.afai` to be masked by `coast_mask.mask`.
"""
function coast_masked!(afai::AFAI, coast_mask::CoastMask)
    afai.afai = coast_masked(afai, coast_mask)
    return nothing
end

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

function plot(afai::AFAI; coast_mask::Union{Nothing, CoastMask} = nothing)
    fig = default_fig()

    if coast_mask === nothing
        afai_data = afai.afai
    else
        afai_data = coast_masked(afai, coast_mask)
    end

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

"""
    afai_median(afai::AFAI)

Return an `Array` of the same size and eltype of `afai.afai` such that the value at each gridpoint
is the median of all `afai.afai` values in a window of size `afai.params.window_size_median_filter`
centered on that gridpoint. `NaN`s are ignored when calculating the median. If every value in the window 
is `NaN`, the value of `afai_median` is also `NaN`.

This can be a lengthy computation.
"""
function afai_median(afai::AFAI)
    lon = afai.lon
    lat = afai.lat
    afai_data = afai.afai
    window_size = afai.params.window_size_median_filter
    afai_median = zeros(eltype(afai_data), size(afai_data))

    for t_i = 1:size(afai_median, 3)

        @info "Compuing AFAI median, time slice $(t_i)/4"

        for lon_i = 1:size(afai_median, 1)
            for lat_i = 1:size(afai_median, 2)
                lons = max(1, lon_i - round(Integer, window_size/2)):min(length(lon), lon_i + round(Integer, window_size/2))
                lats = max(1, lat_i - round(Integer, window_size/2)):min(length(lat), lat_i + round(Integer, window_size/2))
                valid_pixels = [x for x in afai_data[lons, lats, t_i] if !isnan(x)]

                if length(valid_pixels) > 0
                    afai_median[lon_i, lat_i, t_i] = median(valid_pixels)
                else
                    afai_median[lon_i, lat_i, t_i] = NaN
                end
            end
        end
    end

    return afai_median
end

"""
    pixel_classification(afai::AFAI; afai_median::Union{Nothing, Array} = nothing)

Return a `Vector{CartesianIndex{3}}` giving the Sargassum-containing pixels/gridpoints
of `afai.afai`.

This function looks for entries of `afai.afai - afai_median` that are at least as large 
as `afai.params.threshold_median`.

If `afai_median` is not provided, it is computed automatically.
"""
function pixel_classification(afai::AFAI; afai_median::Union{Nothing, Array} = nothing)
    afai_data = afai.afai

    if afai_median === nothing
        afai_median_data = afai_median(afai)
    else
        @assert size(afai_data) == size(afai_median)
        afai_median_data = afai_median
    end

    return findall(x -> x > afai.params.threshold_median, afai_data - afai_median_data)
end

"""
    pixel_unmixing(afai::AFAI; pixel_classification::Union{Nothing, Vector{<:CartesianIndex}} = nothing)

Return an `Array` of the same size and eltype of `afai.afai` such that the value at each gridpoint
is the percentage coverage of Sargassum in that pixel. The coverage is only computed at the pixels 
given by `pixel_classification`, and set to `0.0` elsewhere. 

The coverage is computed as a linear interpolation between global maximum and minimum `afai` values 
provided by `afai.params.afai_U0` and `afai.params.afai_L0`.

If `pixel_classification` is not provided, it is computed automatically.
"""
function pixel_unmixing(afai::AFAI; pixel_classification::Union{Nothing, Vector{<:CartesianIndex}} = nothing)
    afai_data = afai.afai

    if pixel_classification === nothing
        pixel_classification_data = pixel_classification(afai)
    else
        pixel_classification_data = pixel_classification
    end

    unmixed = zeros(eltype(afai_data), size(afai_data))
    afai_U0 = afai.params.afai_U0
    afai_L0 = afai.params.afai_L0

    for pixel_index in pixel_classification_data
        alpha = (afai_data[pixel_index] - afai_L0)/(afai_U0 - afai_L0)
        if alpha > 1
            alpha = 1.0f0
        elseif alpha < 0
            alpha = 0.0f0
        end
        unmixed[pixel_index] = alpha
    end

    return unmixed
end

"""
    coverage(afai::AFAI; unmixed::Union{Nothing, Array} = nothing, step_lon::Integer = 30, step_lat::Integer = 40)


"""
function coverage(afai::AFAI; unmixed::Union{Nothing, Array} = nothing, step_lon::Integer = 30, step_lat::Integer = 40)
    lon = afai.lon
    lat = afai.lat

    if unmixed === nothing
        unmixed_data = pixel_unmixing(afai)
    else
        unmixed_data = unmixed
    end

    lon_bins = Iterators.partition(1:length(lon), step_lon) |> collect
    lat_bins = Iterators.partition(1:length(lat), step_lat) |> collect
    
    coverage_binned = zeros(eltype(unmixed_data), length(lon_bins), length(lat_bins), size(unmixed_data, 3))
    for t = 1:size(coverage_binned, 3)
        for i = 1:length(lon_bins)
            for j = 1:length(lat_bins)
                coverage_binned[i, j, t] = mean(unmixed_data[lon_bins[i], lat_bins[j], t])
            end
        end
    end

    lon_bins_centers = [mean(lon[lon_bin]) for lon_bin in lon_bins]
    lat_bins_centers = [mean(lat[lat_bin]) for lat_bin in lat_bins]
    coverage_tot = sum(coverage_binned[:,:,t] for t = 1:size(coverage_binned, 3))

    return (lon_bins_centers, lat_bins_centers, coverage_tot)
end

"""
    distribution!(coverage::Array; quant::Real = 0.7)
"""
function distribution!(coverage::Array; quant::Real = 0.7)
    thresh = quantile(filter(x -> x > 0.0, coverage), quant)
    for i = 1:size(coverage, 1)
        for j = 1:size(coverage, 2)
            if coverage[i, j] < thresh 
                coverage[i, j] = 0.0
            else
                coverage[i, j] = log10(coverage[i, j])
            end
        end
    end

    coverage = coverage / sum(coverage)

    return nothing
end

