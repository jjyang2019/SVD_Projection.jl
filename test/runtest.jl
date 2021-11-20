
using SVD_Projection
using Test

@testset "SVD_Projection.jl" begin

seed!(1234)
discovery_data = rand([0.,1.,2.],30, 20)
target_data = rand([0.,1.,2.], 30, 15)
(V_discovery, D, V_target) = SVD_Projection.project_left_matrix(discovery_data, target_data, 10)

end
