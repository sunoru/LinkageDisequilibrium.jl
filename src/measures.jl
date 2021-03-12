using LinearAlgebra

_normalize(P, normalized) = normalized ? P : normalize([
    P[1, 1] P[1, 2]
    P[2, 1] P[2, 2]
], 1)

_papb(P) = P[1, 1] + P[1, 2], P[1, 1] + P[2, 1]

_D_max(P, D) = let (pa, pb) = _papb(P)
    D < 0 ?
    max(-pa * pb, -(1 - pa) * (1 - pb)) :
    min(pa * (1 - pb), (1 - pa) * pb)
end

"""
Coefficient of Linkage Disequilibrium

`P` is a matrix for the (not necessary normalized) frequency of haplotypes for two loci A and B.
"""
function D(P; normalized = false)
    P = _normalize(P, normalized)
    pa, pb = _papb(P)
    P[1, 1] - pa * pb
end

"""
A normalized measure of LD.
"""
function D_prime(P; normalized = false)
    P = _normalize(P, normalized)
    vD = D(P; normalized = true)
    vD / _D_max(P, D)
end

"""
The correlation coefficient
"""
function r²(P; normalized = false)
    P = _normalize(P, normalized)
    pa, pb = _papb(P)
    vD = D(P; normalized = true)
    vD ^ 2 / (pa * (1 - pa) * pb * (1 - pb))
end

"""
`snp1` and `snp2` is n-element boolean arrays indicating alleles at two SNPs within n haplotypes.
"""
function r²(snp1::AbstractVector{Bool}, snp2::AbstractVector{Bool})
    P = zeros(2, 2)
    n = length(snp1)
    for i in 1:n
        P[2 - snp1[i], 2 - snp2[i]] += 1
    end
    # P includes numbers of SNP pairs that has following values
    # [(1,1) (1,0)
    #  (0,1) (0,0)]
    r²(P)
end

"""
`snps` is a nxm boolean matrix indicating alleles at m SNPs within n haplotypes.

Return:

`R²`: A matrix of `r²` values for each pair of SNPs.
"""
function r²(snps::Matrix{Bool})
    n, m = size(snps)
    R² = zeros(m, m)
    for i in 1:m-1
        R²[i, i] = 1.0
        for j in i + 1:m
            R²[i, j] = R²[j, i] = r²(snps[:, i], snps[:, j])
        end
    end
    R²
end
