using Dates

"""
    const TREF
"""
const TREF = DateTime(2000, 1)

"""
    time2months
"""
function time2months(yr::Integer, mnth::Integer)
    δyr = yr - year(TREF)
    δmnth = mnth - month(TREF)
    return 12*δyr + δmnth
end

function time2months(yearmonth::Tuple{Integer, Integer})
    δyr = yearmonth[1] - year(TREF)
    δmnth = yearmonth[2] - month(TREF)
    return 12*δyr + δmnth
end

function time2months(time::DateTime)
    δyr = year(time) - year(TREF)
    δmnth = month(time) - month(TREF)
    return 12*δyr + δmnth
end

"""
    months2time
"""
function months2time(months::Integer)
    return yearmonth(TREF + Month(months))
end