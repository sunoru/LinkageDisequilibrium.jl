const testfiles = [
    "ldselect" => "./ldselect.jl",
    "informativeness" => "./informativeness.jl",
]

const args = let x = get(ENV, "TEST_ARGS", "")
    x === "" ? "all" : x
end

const to_test = args === "all" ? [
    (key, value) for (key, value) in testfiles
   ] : [(args, testfiles[findfirst(x -> x[1] === args, testfiles)][2])]

for (key, testfile) in to_test
    @info "Testing $key..."
    @eval module $(gensym())
        include($testfile)
    end
end

