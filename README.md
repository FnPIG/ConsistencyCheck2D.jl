# ConsistencyCheck2D.jl
A Julia script for checking the consistency between a 2D distribution and a 2D possibility distribution using copulas


## Install

Add the required packages

```julia
julia> ]
(v1.0) pkg> add PossibilisticArithmetic
(v1.0) pkg> add BivariateCopulas
(v1.0) pkg> add IntervalArithmetic
(v1.0) pkg> add Test
```
The more up to date version of [PossibilisticArithmetic.jl](https://github.com/AnderGray/PossibilisticArithmetic.jl) may be required:

```julia
julia> ]
(v1.0) pkg> add https://github.com/AnderGray/PossibilisticArithmetic.jl#master
```

## Run

```julia
julia> include("test_minitive_sklar.jl")
```
