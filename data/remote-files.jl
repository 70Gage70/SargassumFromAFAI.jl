using RemoteFiles

function download_data(year::Integer, month::Integer)
    @assert 2017 <= year <+ 2022 "These are the years with full datasets."
    @assert 1 <= month <= 12 "Must be a valid month."

    month_string = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"][month]

    url_start = "https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_atlantic_oceanwatch_AFAI_7D.nc?AFAI"
    url = url_start * "[($(year)-$(month_string)-08T00:00:00Z):7:($(year)-$(month_string)-29T00:00:00Z)][(0.0):1:(38.0)][(-98.0):1:(-38.0)]"

    filename = "afai-$(year)-$(month_string).nc"
    file_symbol = Symbol(filename)
    @RemoteFile(file_symbol, url, file = filename, updates=:never)

    if isfile(file_symbol) && !force
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
