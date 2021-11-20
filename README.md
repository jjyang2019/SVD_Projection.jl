# SVD_Projection.jl

SVD_projection.jl is a Julia module to calculate the left singular matrix of the target genotype data based on the discovery genotype data.

To install the package in Julia v1.0 or later
```{julia}
using Pkg
Pkg.clone("https://github.com/jjyang2019/SVD_projection.jl.git")
```

The function of the modeule is `project_PC` and it requires two data, `A` (discovery) and `B` (target), and one optional argument `N_PC` for the number of principal components calculated (default is 10).

The discovery data and target data are genotype data encoded as 0, 1, or 2 where each row of the represents a variant and each column an subject.  

See the details of using in the Julia help enviornment by typing ?`project_PC`.
Run the test program to learn how the `project_PC` works. This test program  simulates discovery and target data.
```{julia}
using SVD_Projection
include("runtest.jl")
```
