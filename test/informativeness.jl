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

println("Selection with Informativeness Algorithm")

expected = [1, 1, 2, 2, 3]
actual = Int[]

for ξ in 0.2:0.2:1.0

tagSNPs = informative_select(samples; ξ = ξ)
nbins = length(tagSNPs)

println("ξ = $ξ")
println("$nbins tagSNPs: ", tagSNPs)

push!(actual, nbins)

end

@test expected == actual
