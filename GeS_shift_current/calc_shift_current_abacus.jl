using HopTB
using DelimitedFiles
using StaticArrays
using HDF5



function _create_dict_h5(filename::String)
    fid = h5open(filename, "r")
    T = eltype(fid[keys(fid)[1]])
    d_out = Dict{Array{Int64,1}, Array{T, 2}}()
    for key in keys(fid)
        data = read(fid[key])
        nk = map(x -> parse(Int64, convert(String, x)), split(key[2 : length(key) - 1], ','))
        d_out[nk] = permutedims(data)
    end
    close(fid)
    return d_out
end

function create_model_deeph(hamiltonian_path, overlap_path, position_path, lattice_path, orbital_path, spinful)
    orbital_types_f = open(orbital_path, "r")
    orbital_types = Vector{Vector{Int64}}()
    while true
        line_split = split(readline(orbital_types_f))
        if length(line_split) == 0
            break
        end
        orbital_type = parse.(Int64, line_split)
        push!(orbital_types, orbital_type)
    end
    nsites = length(orbital_types)
    site_norbits = (x->sum(x .* 2 .+ 1)).(orbital_types) * (1 + spinful)
    site_norbits_cumsum = cumsum(site_norbits)
    site_norbits_spinless = (x->sum(x .* 2 .+ 1)).(orbital_types)
    site_norbits_cumsum_spinless = cumsum(site_norbits_spinless)
    norbits = sum(site_norbits)
    close(orbital_types_f)

    lattice = readdlm(lattice_path)
    
    hamiltonians = _create_dict_h5(hamiltonian_path)
    overlaps = _create_dict_h5(overlap_path)
    positions = _create_dict_h5(position_path)
    H_type = eltype(hamiltonians[[0, 0, 0, 1, 1]])
    H_R = Dict{SVector{3, Int64}, Matrix{H_type}}()
    S_R = Dict{SVector{3, Int64}, Matrix{Float64}}()
    r_R = Dict{SVector{3, Int64}, SVector{3, Matrix{Float64}}}()
    for key in collect(keys(hamiltonians))
        hamiltonian = hamiltonians[key]
        if (key ∈ keys(overlaps))
            overlap = overlaps[key]
        else
            overlap = zero(hamiltonian)
        end
        if spinful
            overlap = kron([1 0; 0 1], overlap)
        end
        position = Vector{Matrix{Float64}}()
        for direction in 1:3
            key_r = cat(key, direction, dims=1)
            if (key ∈ keys(overlaps))
                position_direction = positions[key_r]
            else
                position_direction = zero(hamiltonian)
            end
            if spinful
                push!(position, kron([1 0; 0 1], position_direction))
            else
                push!(position, position_direction)
            end
        end
        R = key[1:3]; atom_i=key[4]; atom_j=key[5]
        @assert (site_norbits[atom_i], site_norbits[atom_j]) == size(hamiltonian)
        @assert (site_norbits[atom_i], site_norbits[atom_j]) == size(overlap)
        if !(R ∈ keys(H_R))
            H_R[R] = zeros(H_type, norbits, norbits)
            S_R[R] = zeros(Float64, norbits, norbits)
            r_R[R] = SVector{3, Matrix{Float64}}(zeros(Float64, norbits, norbits), zeros(Float64, norbits, norbits), zeros(Float64, norbits, norbits))
        end
        for block_matrix_i in 1:site_norbits_spinless[atom_i]
            for block_matrix_j in 1:site_norbits_spinless[atom_j]
                if spinful
                    index_i = site_norbits_cumsum_spinless[atom_i] - site_norbits_spinless[atom_i] + block_matrix_i
                    index_j = site_norbits_cumsum_spinless[atom_j] - site_norbits_spinless[atom_j] + block_matrix_j
                    H_R[R][index_i, index_j] = hamiltonian[block_matrix_i, block_matrix_j]
                    H_R[R][index_i + norbits ÷ 2, index_j] = hamiltonian[block_matrix_i + site_norbits_spinless[atom_i], block_matrix_j]
                    H_R[R][index_i, index_j + norbits ÷ 2] = hamiltonian[block_matrix_i, block_matrix_j + site_norbits_spinless[atom_j]]
                    H_R[R][index_i + norbits ÷ 2, index_j + norbits ÷ 2] = hamiltonian[block_matrix_i + site_norbits_spinless[atom_i], block_matrix_j + site_norbits_spinless[atom_j]]

                    S_R[R][index_i, index_j] = overlap[block_matrix_i, block_matrix_j]
                    S_R[R][index_i + norbits ÷ 2, index_j] = overlap[block_matrix_i + site_norbits_spinless[atom_i], block_matrix_j]
                    S_R[R][index_i, index_j + norbits ÷ 2] = overlap[block_matrix_i, block_matrix_j + site_norbits_spinless[atom_j]]
                    S_R[R][index_i + norbits ÷ 2, index_j + norbits ÷ 2] = overlap[block_matrix_i + site_norbits_spinless[atom_i], block_matrix_j + site_norbits_spinless[atom_j]]

                    for direction in 1:3
                        r_R[R][direction][index_i, index_j] = position[direction][block_matrix_i, block_matrix_j]
                        r_R[R][direction][index_i + norbits ÷ 2, index_j] = position[direction][block_matrix_i + site_norbits_spinless[atom_i], block_matrix_j]
                        r_R[R][direction][index_i, index_j + norbits ÷ 2] = position[direction][block_matrix_i, block_matrix_j + site_norbits_spinless[atom_j]]
                        r_R[R][direction][index_i + norbits ÷ 2, index_j + norbits ÷ 2] = position[direction][block_matrix_i + site_norbits_spinless[atom_i], block_matrix_j + site_norbits_spinless[atom_j]]
                    end
                else
                    index_i = site_norbits_cumsum[atom_i] - site_norbits[atom_i] + block_matrix_i
                    index_j = site_norbits_cumsum[atom_j] - site_norbits[atom_j] + block_matrix_j
                    H_R[R][index_i, index_j] = hamiltonian[block_matrix_i, block_matrix_j]
                    S_R[R][index_i, index_j] = overlap[block_matrix_i, block_matrix_j]
                    for direction in 1:3
                        r_R[R][direction][index_i, index_j] = position[direction][block_matrix_i, block_matrix_j]
                    end
                end
            end
        end
    end

    tb_model = TBModel{H_type}(norbits, lattice, isorthogonal=false)
    tb_model.site_norbits = site_norbits
    tb_model.nsites = nsites
    tb_model.isspinful = spinful
    
    tb_model.hoppings = H_R
    tb_model.overlaps = S_R
    tb_model.positions = r_R

    set_orbital_types!(tb_model, orbital_types, isspinful=spinful)
    return tb_model
end


fermi_level = 0.3970374823289725

tb_model_from_DeepH = create_model_deeph(
"./h5files/hamiltonians.h5",
"./h5files/overlaps.h5",
"./h5files/positions.h5",
"./h5files/lat.dat",
"./h5files/orbital_types.dat",
false
)

tb_model = HopTB.SharedTBModel(tb_model_from_DeepH)

alpha = 1
beta = 1
gamma = 1
ωs = collect(range(0, stop=5, length=500))
kmesh = [200, 200, 1]
gauss_width = 0.02
sc = HopTB.Optics.get_shift_cond(tb_model, alpha, beta, gamma, ωs, fermi_level, kmesh, ϵ=gauss_width)
open(string("shiftcond_", alpha, beta, gamma, ".dat"), "w") do f
    writedlm(f, [ωs sc])
end
