module SVD_Projection

export project_left_matrix
using LinearAlgebra

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
    return U,D,V
end

function cal_snp_freq(obj)
    res = zeros(eltype(obj), size(obj,1))
    for i in 1:size(obj,1)
        res[i] = cal_p(obj[i,:])
    end
    return res
end

function cal_V_projection(obj, U, D)
    transpose(obj) * U * inv(Diagonal(D))
end

function project_left_matrix(train, test, N_PC=10)
    #println("Calculate SNP freq")
    snp_freq = cal_snp_freq(train)

    #println("Normalization Train data")
    raw_norm = normalized_data(train, snp_freq);

    #println("Calculate SVD")
    (U,D,V_train) = my_svd(raw_norm, N_PC);

    #println("Normalization Test Data")
    test_norm = normalized_data(test , snp_freq)

    #println("Project Test Data to V matrix")
    V_test = cal_V_projection(test_norm, U, D)
    return (V_train, D, V_test)
end


end
