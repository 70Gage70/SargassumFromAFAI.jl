using RemoteFiles

"""
    download_data(year::Integer, month::Integer)

Download the AFAI data for the given year and month from 

`https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_atlantic_oceanwatch_AFAI_7D.html`

Can be applied as `download_data(year::Integer, month::Vector{<:Integer})` to download multiple months in a given 
year at once.

Four files are downloaded, one on the 8th, 15th, 22nd and 29th of each month, each of which are 7-day aggregates. 

The files are named "afai-year-month.nc" and are stored in the `data` folder of this package. 

Use the function `data_path(year, month)` to obtain the path to the file. And `data_rm(year, month)` to remove it.
"""
function download_data(year::Integer, month::Integer)
    @assert 2017 <= year <+ 2022 "These are the years with full datasets."
    @assert 1 <= month <= 12 "Must be a valid month."

    month_string = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"][month]

    url_start = "https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_atlantic_oceanwatch_AFAI_7D.nc?AFAI"
    url = url_start * "[($(year)-$(month_string)-08T00:00:00Z):7:($(year)-$(month_string)-29T00:00:00Z)][(0.0):1:(38.0)][(-98.0):1:(-38.0)]"

    filename = "afai-$(year)-$(month_string).nc"
    file_symbol = Symbol(filename)
    @RemoteFile(file_symbol, url, file = filename, updates=:never)

    if isfile(file_symbol) # && !force
        @info "Already downloaded at $(path(file_symbol))."
        return nothing
    end

    RemoteFiles.download(RemoteFiles.Http(), file_symbol)

    return nothing
end

function download_data(year::Integer, months::Union{Vector{<:Integer}, AbstractRange})
    for month in months
        download_data(year, month)
    end

    return nothing
end

"""
    data_path(year::Integer, month::Integer)

Return the path to the raw data set for the given month and year, if it exists.
"""
function data_path(year::Integer, month::Integer)
    @assert 2017 <= year <+ 2022 "These are the years with full datasets."
    @assert 1 <= month <= 12 "Must be a valid month."

    month_string = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"][month]

    url_start = "https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_atlantic_oceanwatch_AFAI_7D.nc?AFAI"
    url = url_start * "[($(year)-$(month_string)-08T00:00:00Z):7:($(year)-$(month_string)-29T00:00:00Z)][(0.0):1:(38.0)][(-98.0):1:(-38.0)]"

    filename = "afai-$(year)-$(month_string).nc"
    file_symbol = Symbol(filename)
    @RemoteFile(file_symbol, url, file = filename, updates=:never)

    if !(isfile(file_symbol))
        @info "File has not been downloaded."

        return nothing
    end

    return path(file_symbol)
end

"""
    data_path(year::Integer, month::Integer)

Remove the raw data set for the given month and year, if it exists.
"""
function data_rm(year::Integer, month::Integer)
    @assert 2017 <= year <+ 2022 "These are the years with full datasets."
    @assert 1 <= month <= 12 "Must be a valid month."

    month_string = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"][month]

    url_start = "https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_atlantic_oceanwatch_AFAI_7D.nc?AFAI"
    url = url_start * "[($(year)-$(month_string)-08T00:00:00Z):7:($(year)-$(month_string)-29T00:00:00Z)][(0.0):1:(38.0)][(-98.0):1:(-38.0)]"

    filename = "afai-$(year)-$(month_string).nc"
    file_symbol = Symbol(filename)
    @RemoteFile(file_symbol, url, file = filename, updates=:never)

    if !(isfile(file_symbol))
        @info "File has not been downloaded."

        return nothing
    end

    rm(file_symbol)

    return nothing
end
