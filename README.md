# L2OViz.jl
L2OViz.jl visualizes the solutions to multiple instances of an optimization problem.
It supports visualizing the solutions of the same instances from multiple solvers for comparison.
It also has a special feature for visualizing (the most interesting rows and columns of) matrix variables in a 2D layout of subplots which correspond to coordinates in the matrix.
This can be useful when the variables correspond to a graph for example.



## Data Format Specifications
For each variable, the values should be stored in a `Matrix` where each row contains the values of the variable in each problem instance.
For example,
```Julia
y_dim = 5
n_instances = 3
y = randn(y_dim, n_instances)
```
contains the values of a 5-dimensional variable for 3 instances of an optimization problem.

`plot_variable` accepts a variable number of `Matrix` inputs, each corresponding to a solver.

### Matrix variables
Matrix variable data should be provided in COO format.
It is assumed that, across all the problem instances, the same matrix variable has the same dimensions and sparsity structure.
Therefore, the values are still stored as a `Matrix` with `n_instances` columns, and each column contains the nonzero values of the matrix variable in each problem instance.
In other words, all the variables are treated as vector variables in L2OViz.jl.

`plot_matrix_variable` accepts a variable number of `Matrix` inputs, each corresponding to a solver.


## Visualization
The values of each variable entry across all the problem instances are visualized in a scatter point subplot.
`plot_variable` simply places the subplots side-by-side.
`plot_matrix_variable` arranges the subplots into a grid layout, where the subplot at coordinate `(i, j)` visualizes the `(i, j)` entry of the variable as specified in `(I, J)`.

The data of Solver A and Solver B do not have to be for the same problem instances.
In this case, different `x` should be provided.


### Thresholding
When the dimension of the variable to visualize is too high, `vis_threshold` limits the number of entries that are visualized.
`significance_fn` is used to select the most interesting entries of the variable.


## Example: Optimal Power Flow
`exp/viz_opf.jl` is a utility script for visualizing OPF solution data stored in HDF5 files, using system data from PGLib.jl and PowerModels.jl.

**Arguments**

| Flag | Description | Default |
|---|---|---|
| `--system` | PGLib system name (e.g. `pglib_opf_case14_ieee`) | *(required)* |
| `--h5` | Path to the `.h5` solution file for each solver | *(required)* |
| `--variables` | Variable names to visualize; one image is saved per variable | *(required)* |
| `--solver-names` | Display name for each solver, matching the order of `--h5` | `Solver 1`, `Solver 2`, … |
| `--output-dir` | Directory where output images are saved | `.` |
| `--vis-threshold` | Max entries / rows+columns to show per variable | `20` |
| `--flat` | Force `plot_variable` for all variables | off |

**Variable dispatch**

By default (without `--flat`), the script infers the plot type from the variable dimension:
- Dimension equals the number of **branches** → `plot_matrix_variable`, with `I` and `J` taken from `f_bus` and `t_bus` of the network branches in sorted key order.
- Dimension equals the number of **buses** → `plot_variable`.

With `--flat`, all variables use `plot_variable` regardless of dimension, and the network is not loaded.

**x-axis**

Currently, the x-axis is `sum(pd .+ qd, dims=1)` (the total active and reactive load per problem instance).

**Usage**
```bash
julia --project=exp exp/viz_opf.jl \
    --system pglib_opf_case14_ieee \
    --h5 solver_a.h5 solver_b.h5 \
    --solver-names "Solver A" "Solver B" \
    --variables pg vm \
    --output-dir results/
```

Each `.h5` file must contain a dataset per variable, and `pd` and `qd` datasets.
Output images are named `{system}_{variable}.png`.

