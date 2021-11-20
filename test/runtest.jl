using SVD_Projection
using Test, Random
Random.seed!(1234)
discovery = rand([0., 1., 2.], 30, 20);
target= rand([0., 1., 2.], 30, 15);

V_PC = project_PC(discovery, target, N_PC=3)
@test V_PC[1,:] == [-0.12869065723376472; -0.05776618956235076; -0.12278055691277556]
