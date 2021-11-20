# SVD_Projection.jl

SVD_projection.jl is a Julia module to calculate the left singular matrix of the target genotype data based on the discovery genotype data.

To install the package in Julia v1.0 or later
```{julia}
using Pkg
Pkg.clone("https://github.com/jjyang2019/SVD_projection.jl.git")
```

The training(discovery), testing(target) data, and test_project_V_matrix.jl are in the test folder to demonstrate the usage. Each row of the training data represents a variant and each column an subject. Each cell in the training data matrix is coded as the number of the effect allele. The test data has the same format. To run the test program, simply type
```{julia}
using SVD_Projection
include("runtest.jl")
```
