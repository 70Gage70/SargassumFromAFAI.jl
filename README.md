# SargassumFromAFAI

This is a [Julia](https://julialang.org/) package for computing *Sargassum* coverage distributions from AFAI (Alternate Floating Algae Index). The  data are generated based on the [7-day cumulative USF AFAI Fields](https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_atlantic_oceanwatch_AFAI_7D.html) and according to the analysis presented in [Wang and Hu, 2016](https://www.sciencedirect.com/science/article/abs/pii/S0034425716301833). It should be noted  that the analysis performed here is *not* identical to Wang and Hu, 2016 and the results here are *not* identical to the results shown in the  [*Sargassum* Watch System](https://optics.marine.usf.edu/projects/saws.html).

These data and results are presented for scientific research purposes **only**.

## Installation

In the Julia REPL, execute the following

```julia
import Pkg
Pkg.add("https://github.com/70Gage70/SargassumFromAFAI.jl.git")
```

## Quickest Start

```julia

using SargassumFromAFAI


```


## Quickest Start