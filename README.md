# GNRProdEst

This package implements the production function estimaton method from Gandhi, Navarro, Rievers (2020).

## Installation

Use Julia's integrated package manager [Pkg.jl](https://github.com/JuliaLang/Pkg.jl)


```Julia
using Pkg
Pkg.add("GNRProdEst")
```

## Usage

You can use estimation routine in two ways. Either run both stages using `gnrprodest` or run each estimation stage separately with `gnrfirststage` and `gnrsecondstage`. Versions of these functions that end with an `!` modify their inputs (most importantly the dataframe).
### Arguments
All arguments in any of the functions are keyword arguments. The necessary arguments are:
- `data` is a DataFrame that includes all columns defined by other arguments.
- `output` is the Symbol and represents the firms natural logarithm of gross output.
- `flexible_input` is a Symbol representing the flexible input of the production function.
- `fixed_input` is a Symbol or a vector of Symbols representing any fixed production input.
- ln_share_flex_y



                  flexible_input::Symbol, 
                  fixed_inputs::Union{Symbol,Array{Symbol}}, 
                  ln_share_flex_y_var::Symbol = :NotDefinedByUser, 
                  id::Symbol, 
                  time::Symbol, 
                  fes_starting_values::Vector = [missing],
                  ses_starting_values::Vector = [missing],
                  share_degree::Int = 3,
                  lm_tfp_degree::Int = 3,
                  int_const_series_degree::Int = 3,
                  opts::Dict = Dict())


```Julia
using GNRProdEst

gnrprodest(data = )
```

## Example


## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

[MIT](https://choosealicense.com/licenses/mit/)