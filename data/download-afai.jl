using Downloads

year = 2018
months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]

url_start = "https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_atlantic_oceanwatch_AFAI_7D.nc?AFAI"
urls = String[]

for i = 1:length(months)
    push!(urls, url_start * "[($(year)-$(months[i])-08T00:00:00Z):7:($(year)-$(months[i])-29T00:00:00Z)][(0.0):1:(38.0)][(-98.0):1:(-38.0)]")
end

for i = 1:length(months)
    Downloads.download(urls[i], "afai-$(year)-$(months[i]).nc")
end

