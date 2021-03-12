module LinkageDisequilibrium

include("./utils.jl")

include("./measures.jl")

export fisher_exact_test, χ²_test
include("./hypothesis_tests.jl")

export ld_select
include("./ld_select.jl")

export informative_select
include("./informative_select.jl")

end
