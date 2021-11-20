module SVD_Projection
using LinearAlgebra
export project_PC

function cal_p(x)
    sum = 0.
    num = 0
    for z in x
        if ! ismissing(z)
            num += 1
            if z == 2
                sum +=1
            elseif z == 1
                sum += 0.5
            end
        end
    end
    res = sum / num
    return res
end

function valid_num_ind(x)
    isfinite.(x) .& .!ismissing.(x) .& .!isnan.(x)
end

function normalized_data(obj,snp_freq)
    res = zeros(eltype(obj),size(obj))
    for i in 1:size(obj,1)
        x = obj[i,:]
        p = snp_freq[i]
        z = (x .- 2 * p) ./ sqrt(2 * p * (1.0 - p))
        z[.!valid_num_ind(z)] .= 0.
        res[i,:] = z
    end
    res = res ./ sqrt(size(res,1))
    return res
end

function get_V_matrix(obj, pc_num)
    myfit = eigen(transpose(obj) * obj, sortby = x -> -abs(x))

    if pc_num < size(obj,2)
        eigenval = myfit.values[1:pc_num]
        eigenvec = myfit.vectors[:,1:pc_num]
    else
        eigenval = myfit.values
        eigenvec = myfit.vectors
    end
    return eigenval, eigenvec
end

function my_svd(obj, pc_num)
    if pc_num > size(obj, 2)
        pc_num = size(obj, 2)
    end
    res = get_V_matrix(obj, pc_num)

    D = sqrt.(res[1])
    V = res[2]

    U = obj * V * inv(Diagonal(D))
    return U, D, V
end

function cal_snp_freq(obj)
    res = zeros(eltype(obj), size(obj,1))
    for i in 1:size(obj,1)
        res[i] = cal_p(obj[i,:])
    end
    return res
end

function cal_V_projection(obj, U, D)
    res = transpose(obj) * U * inv(Diagonal(D))
    return res
end

"""
    project_PC(A, B; N_PC=10)

Calculate the PCs (the right projected singular matrix) of `B`.
Based on the singular value decomposition, the matrix `A` is represented as `A=USV'`.
Therefore, the PCs are columns of matrix `V`. To project `B` onto `V`, the `project_PC`
return `B'US'`.

# Arguments
* `A`: The discovery genotype matrix, where eeach column is a subject and each row is a variant.
* `B`: The target genotype matrix, where eeach column is a subject and each row is a variant.
* `N_PC`: The optional argument to determine the number of columnes in the returned PC matrix.

# Examples
```julia
julia> using Random
julia> Random.seed!(1234)
julia> discovery = rand([0.,1.,2.], 30, 20);
julia> target = rand([0.,1.,2.], 30, 15);
julia> target_PC = project_PC(discovery, target, N_PC=3)
```
"""
function project_PC(A, B; N_PC=10)
    snp_freq = cal_snp_freq(A)
    raw_norm = normalized_data(A, snp_freq);
    (U, D, V) = my_svd(raw_norm, N_PC);
    test_norm = normalized_data(B, snp_freq)
    res_PC = cal_V_projection(test_norm, U, D)
    return res_PC
end
end
