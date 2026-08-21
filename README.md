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
### Using `gnrprodest` and `gnrprodest!`
Both versions have exactly the same arguments. The difference is, that gnrprodest modifies its arguments (especially `data`). `fes_returns` is a NamedTuple of first stage return objects, `ses_returns` is the equivalent for the second stage, `results_table` is a DataFrame printed at the end containing key statistics, and `estimation_sample` is the DataFrame used in the estimation routine.
```Julia
using GNRProdEst

fes_returns, ses_returns, results_table, estimation_sample = gnrprodest(data = DataFrame, output = :revenue, flexible_input = :flexible_input, fixed_inputs = [:fixed_input1 :fixed_input2 ...], id = :id, time = :time)

fes_returns, ses_returns, results_table = gnrprodest!(data = DataFrame, output = :gross_output, flexible_input = :flexible_input, fixed_inputs = [:fixed_input1 :fixed_input2 ...], id = :id, time = :time)
```
#### Arguments
All arguments in any of the functions are keyword arguments. The necessary arguments are:
- `data` is a DataFrame that includes all columns defined by other arguments.
- `output` is the Symbol and represents the firms natural logarithm of gross output.
- `flexible_input` is a Symbol representing the natural logarithm of the flexible input of the production function.
- `fixed_inputs` is a Symbol or a vector of Symbols representing any natural logarithm of fixed production inputs.
- `id` is a Symbol of each firm's unique identifier.
- `time` is a Symbol identifying the observation`s time period.

#### Optional Arguments
- `ln_share_flex_y` is a Symbol for the natural logarithm of the flexible input share of revenue. If it is not defined by the user, the package calculates it internally.
- `share_degree` is the degree of the polynomial series in the share regression. Its default value is `3` if it is not defined by the user.
- `lm_tfp_degree` is the degree of the law of motion for $\omega$ (tfp) over time. The default value if not defined by user is `3`.
- `int_const_series_degree` is the degree of the polynomial series of the constant of integration in the second stage. Its default value is `3` if not defined by the user.
- `fes_starting_values` is a vector of starting values for the first stage non-linear least squares regression. If the user does not define them, the package chooses them using a linear regression.
- `ses_starting_values` is a vector of starting values for the second stage GMM estimation routine. If the user does not define them, the package chooses them using a linear regression.
- `boot_reps` is an integer indicating the number of repetitions the Bootstrap algorithm uses in the standard error estimation. The default is 200.
  Bootstrapping is by far the slowest part of the estimation because every repetition re-estimates both stages. The repetitions run on Julia's threads, so starting Julia with several threads speeds this up considerably (e.g. `julia --threads=auto`, or set the `JULIA_NUM_THREADS` environment variable). With a single thread they run one after another.
- `opts` is a dictionary of futher options.
    - `fes_print_starting_values` prints first stage starting values if set to `true`.
    - `fes_print_results` prints first stage results if set to `true`.
    - `fes_optimizer` can be used to change the numerical optimization algorithm in the first stage. The default is NelderMead. You can find a list of available optimizers with an in-dept description here: [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/)
    - `fes_optimizer_options` allows the user to configure the options of the first stage numerical optimization routine. Use an `Optim.Options()` object to define them. The configurable options depend on the optimizer chosen. You can find a list of options here: [Optim.Options](https://julianlsolvers.github.io/Optim.jl/stable/user/config/)
    - `ses_print_starting_values` prints first stage starting values if set to `true`.
    - `ses_print_results` prints first stage results if set to `true`.
    - `ses_optimizer` can be used to change the numerical optimization algorithm in the second stage. The default is NelderMead. You can find a list of available optimizers with an in-dept description here: [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/)
    - `ses_optimizer_options` allows the user to change the options of the second stage numerical optimization routine. Use an `Optim.Options()` object to define them. The configurable options depend on the optimizer chosen. You can find a list of options here: [Optim.Options](https://julianlsolvers.github.io/Optim.jl/stable/user/config/)
    - `print_results` prints a summary table of key statistics including (among others) standard errors and confidence intervals. Deflaut is `true`.

### Using `gnrfirststage` and `gnrfirststage!`
Instead of using one command for both stages, you can estimate both stages independently. The command for the estimation of the first stage is `gnrfirststage` or `gnrfirststage!`. The difference is that the latter changes its inputs (especially the data). See the list of argument descriptions above.`fes_returns` is a NamedTuple of first stage return objects and `estimation_sample` is the DataFrame used in the estimation routine.

```Julia
using GNRProdEst

fes_returns, estimation_sample = gnrfirststage(data = DataFrame, output = gross_output, flexible_input = :flexible_input, fixed_inputs = [:fixed_input1 :fixed_input2 ...],ln_share_flex_y = :ln_share_flex_y, share_degree = :share_degree, starting_values = :starting_values [, opts = Dict])

fes_returns = gnrfirststage!(data = DataFrame, output = gross_output, flexible_input = :flexible_input, fixed_inputs = [:fixed_input1 :fixed_input2 ...],ln_share_flex_y = :ln_share_flex_y, share_degree = :share_degree, starting_values = :starting_values [, opts = Dict])
```

### Using `gnrsecondstage` and `gnrsecondstage!`
The command for the estimation of the second stage is `gnrsecondstage` or `gnrsecondstage!`. The difference is that the latter changes its inputs (especially the data). See the list of argument descriptions above. The only additional argument is `fes_returns` representing the object returned from the first stage estimation. Use `opts` as described above to modify additional options, such as the numerical optimization algorithm. `ses_returns` is a NamedTuple of return objects from the second stage estimation.

```Julia
using GNRProdEst

ses_returns, estimation_sample = gnrsecondstage(data = DataFrame,flexible_input = :flexible_input, fixed_inputs = [:fixed_input1 :fixed_input2 ...], id = :id, time = :time, fes_returns = fes_returns, int_const_series_degree = int_const_series_degree, lm_tfp_degree = lm_tfp_degree, starting_values = ses_starting_values [, opts = Dict])

ses_returns = gnrsecondstage!(data = DataFrame,flexible_input = :flexible_input, fixed_inputs = [:fixed_input1 :fixed_input2 ...], id = :id, time = :time, fes_returns = fes_returns, int_const_series_degree = int_const_series_degree, lm_tfp_degree = lm_tfp_degree, starting_values = ses_starting_values [, opts = Dict])
```

## Example
This example uses 500 firms from [Gandhi, Navarro, Rivers (2020)](https://www.journals.uchicago.edu/doi/abs/10.1086/707736?af=R&mobileUi=0&)'s replication package data. It estimates a gross output production function with one flexible input $i$ and one fixed input $k$.
```Julia

# Load required dependencies (use 'using Pkg' and Pkg.add("DataFrames"), Pkg.add("GNRProdEst") or Pkg.add("CSV") if you do not have any of these packages in your environment.)
using GNRProdEst, DataFrames, CSV

# Read in replication data
rep_data = CSV.read("test/GNR_data_500.csv", DataFrame)

# Run both estimation stages at the same time
gnr_res = GNRProdEst.gnrprodest!(data = rep_data, 
                                output = :yg, 
                                flexible_input = :i, 
                                fixed_inputs = :k, 
                                ln_share_flex_y = :si, 
                                id = :id, 
                                time = :time);


```

Output should be
```Julia
.......................................................................................
.......................................................................................
..........................
Number of observations: 14500
┌────────────────────┬──────────┬─────────┬───────────┬─────────┬──────────┬───────────┐
│           Variable │     Stat │      SE │    t_stat │ P_Value │ Conf_low │ Conf_high │
├────────────────────┼──────────┼─────────┼───────────┼─────────┼──────────┼───────────┤
│ mean elasticity: k │  0.28670 │ 0.04377 │   6.54973 │ 0.00000 │  0.20091 │   0.37250 │
│ mean elasticity: i │  0.62627 │ 0.00151 │ 414.77249 │ 0.00000 │  0.62331 │   0.62923 │
│               cons │  0.16803 │ 0.03493 │   4.81045 │ 0.00000 │  0.09956 │   0.23649 │
│                  ω │  0.76907 │ 0.05542 │  13.87800 │ 0.00000 │  0.66045 │   0.87768 │
│                 ω² │  0.06544 │ 0.06752 │   0.96927 │ 0.33241 │ -0.06689 │   0.19778 │
│                 ω³ │ -0.03981 │ 0.02510 │  -1.58593 │ 0.11276 │ -0.08900 │   0.00939 │
└────────────────────┴──────────┴─────────┴───────────┴─────────┴──────────┴───────────┘
```

Alternatively, you can estimate both stages sparately with
```Julia
# Load required dependencies (use 'using Pkg' and Pkg.add("DataFrames"), Pkg.add("GNRProdEst") or Pkg.add("CSV") if you do not have any of these packages in your environment.)
using GNRProdEst, DataFrames, CSV

# Read in replication data
rep_data = CSV.read("test/GNR_data_500.csv", DataFrame)

# Define some options to print results
opts = Dict("fes_print_results" => true, "ses_print_results" => true)

# Run the first estimation stage
fes_return_obj, est_sample = GNRProdEst.gnrfirststage(data = rep_data,
                                                      output = :yg,
                                                      flexible_input = :i,
                                                      fixed_inputs = :k,
                                                      ln_share_flex_y = :si,
                                                      opts = opts);

```

Output should be

```Julia
First stage results:
┌──────────┬──────────┬──────────┐
│ Variable │        γ │       γ' │
├──────────┼──────────┼──────────┤
│ constant │  0.65239 │  0.67590 │
│        k │ -0.00146 │ -0.00152 │
│      k⋅k │ -0.00050 │ -0.00052 │
│    k⋅k⋅k │ -0.00033 │ -0.00035 │
│    k⋅k⋅i │  0.00111 │  0.00115 │
│      k⋅i │  0.00123 │  0.00128 │
│    k⋅i⋅i │ -0.00112 │ -0.00116 │
│        i │ -0.02119 │ -0.02195 │
│      i⋅i │  0.00476 │  0.00494 │
│    i⋅i⋅i │ -0.00005 │ -0.00005 │
└──────────┴──────────┴──────────┘
```

Afterwards, you can run the second estimation step with
```Julia

ses_res = GNRProdEst.gnrsecondstage(data = est_sample,
                                    flexible_input = :i,
                                    fixed_inputs = :k,
                                    id = :id,
                                    time = :time,
                                    fes_returns = fes_return_obj,
                                    opts = opts);
```
The output should be
```Julia
Integration constant series parameters
┌──────────┬──────────┐
│ Variable │ Estimate │
├──────────┼──────────┤
│        k │  0.38825 │
│      k⋅k │ -0.02485 │
│    k⋅k⋅k │  0.00247 │
└──────────┴──────────┘
All output elasticities:
┌──────────┬─────────┬─────────┬─────────┬─────────┐
│ Variable │    Mean │      SD │     Min │     Max │
├──────────┼─────────┼─────────┼─────────┼─────────┤
│        k │ 0.28670 │ 0.01767 │ 0.27269 │ 0.38782 │
│        i │ 0.62627 │ 0.00454 │ 0.60449 │ 0.65235 │
└──────────┴─────────┴─────────┴─────────┴─────────┘
Productivity:
┌──────────┬─────────┬─────────┬──────────┬─────────┐
│ Variable │    Mean │      SD │      Min │     Max │
├──────────┼─────────┼─────────┼──────────┼─────────┤
│        ω │ 0.80983 │ 0.34225 │ -0.47409 │ 1.67487 │
│        Ω │ 2.38171 │ 0.82280 │  0.62245 │ 5.33812 │
│        v │ 0.80987 │ 0.43333 │ -0.86880 │ 2.30456 │
└──────────┴─────────┴─────────┴──────────┴─────────┘
Productivity Law of Motion:
┌──────────┬──────────┐
│ Variable │ Estimate │
├──────────┼──────────┤
│ constant │  0.16803 │
│        ω │  0.76907 │
│       ω² │  0.06544 │
│       ω³ │ -0.03981 │
└──────────┴──────────┘
```
## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
Please cite this package if you use it for published research.

[MIT](https://choosealicense.com/licenses/mit/)
