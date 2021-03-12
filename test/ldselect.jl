using Test

using LinkageDisequilibrium

samples = Bool[
    1 1 1 1 1 1 1 1
    0 1 1 1 0 1 0 1
    0 1 0 1 0 1 0 1
    1 1 0 1 0 1 0 0
    1 1 0 1 1 1 0 1
    0 0 1 1 1 1 0 0
    0 0 0 0 0 0 0 0
]

println("LD-Select Algorithm")

expected = [2, 5, 5, 6, 7, 7, 7]
actual = Int[]

for τ in 0.2:0.1:0.8

bins = ld_select(samples; r²_threshold = τ)
nbins = length(bins)
tagSNPs = collect(keys(bins))

println("τ = $τ")
println("$nbins bins: ", tagSNPs)

push!(actual, nbins)

end

@test expected == actual
