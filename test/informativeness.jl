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

expected = [1, 1, 2, 2, 3]
actual = [
    length(informative_select(samples; ξ = ξ))
    for ξ in 0.2:0.2:1.0
]

@test expected == actual
