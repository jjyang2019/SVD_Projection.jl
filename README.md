# SVD_Projection.jl

SVD_projection.jl is a Julia module to calculate the projected principal components (PCs) of the target genotype data based on the discovery genotype data.

To install the package in Julia v1.0 or later
```{julia}
using Pkg
Pkg.add(url="https://github.com/jjyang2019/SVD_projection.jl.git")
```

There is only one function in  this modeule which is `project_PC(A, B, N_PC)`. The two required arguments are `A` (discovery) and `B` (target), and one optional argument `N_PC`.

The discovery data (`A`) and target data (`B`) are genotype data encoded as 0, 1, or 2 where each row  represents a variant and each column represents a subject. To convert plink bfile to discovery or target data, one can use the plink command `plink --bfile bfile --export A-transpose --out out_file`.  

See the details of using `project_PC` in the Julia help enviornment by typing ?`project_PC`.
Run the test program to learn how the `project_PC` works.  
```{julia}
using SVD_Projection
include("runtest.jl")
```
