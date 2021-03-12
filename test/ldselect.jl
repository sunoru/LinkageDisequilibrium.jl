using Test

using LinkageDisequilibrium

samples = Bool[
    0 1 1 1 0 1 0 1
    0 1 0 1 0 1 0 1
    1 1 0 1 0 1 0 0
    1 1 0 1 1 1 0 1
    1 1 1 1 1 1 1 1
    0 0 1 1 1 1 0 0
    0 0 0 0 0 0 0 0
]

expected = [2, 5, 5, 6, 7, 7, 7]
actual = [
    length(ld_select(samples; r²_threshold = τ))
    for τ in 0.2:0.1:0.8
]

@test expected == actual
