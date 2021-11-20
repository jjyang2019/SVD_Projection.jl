
using SVD_Projection
using Test
using DataFrames, CSV

y = myf(2.0)
println(y)
train = rand([0,1,2],30, 20)
test = rand([0,1,2], 30, 15)
(V_train, D, V_test) = projectV(train, test, 10)

@testset "SVD_Projection.jl" begin

function read_geno(f; header=false)
    df = CSV.read(f, DataFrame, header=header)
    res = Matrix{Union{Missing, Float64}}(df)
    return res
end

IN_train="train_data.csv"
IN_test="test_data.csv"

println("Read Train Data: ", IN_train)
train = read_geno(IN_train)

println("Read Test Data: ", IN_test)
test_data = read_geno(IN_test)

num_pc = 12 # number of principal components extracted
(V_train, D, V_test) = project1(train, test, num_pc)


println("Write output")
CSV.write("V_train.csv", DataFrame(V_train), header=false)
CSV.write("D_train.csv", DataFrame([D]), header=false)
CSV.write("V_test.csv", DataFrame(V_test), header=false)

end
