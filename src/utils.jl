# Assuming all SNPs are biallelic.
function make_binary(haplotypes::AbstractVector{T}) where T <: AbstractVector
    n = length(haplotypes)
    n ≡ 0 && return zeros(Bool, 0, 0)
    h1 = haplotypes[1]
    m = length(h1)
    snps = zeros(bool, n, m)
    m ≡ 0 && return snps
    for i in 2:n
        h = haplotypes[i]
        for j in 1:m
            snps[i, j] = h[j] ≢ h1[j]
        end
    end
    snps
end

# Calculate frequencies and make major alleles `true`s.
function prepare_select(haplotypes::Matrix{Bool})
    n, m = size(haplotypes)
    haplotypes = copy(haplotypes)
    fs = vec(sum(haplotypes; dims = 1) / n)
    for j in 1:m
        if fs[j] < 0.5
            haplotypes[:, j] .= .!haplotypes[:, j]
            fs[j] = 1 - fs[j]
        end
    end
    fs, haplotypes
end
