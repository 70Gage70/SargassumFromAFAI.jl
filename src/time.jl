"""
    const TREF

The reference time for the `time` variable in distribution NetCDFs. 

Equal to `January 1, 2000`.
"""
const TREF = DateTime(2000, 1)

"""
    time2months(yr::Integer, mnth::Integer)

Calculate the number of months since [`TREF`](@ref) of the date with year `yr` and month `mnth`.

    time2months(yearmonth::Tuple{Integer, Integer})

Compute `time2months(yearmonth[1], yearmonth[2])`.

    time2months(time::DateTime)

Compute `time2months(year(time), month(time))`.

### Example

```julia
time2months(2018, 4)
219
```
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
    months2time(months::Integer)

Calculate the date in `(year, month)` format after `months` months have elapsed since [`TREF`](@ref).

### Example

```julia
months2time(219)
(2018, 4)
```
"""
function months2time(months::Integer)
    return yearmonth(TREF + Month(months))
end