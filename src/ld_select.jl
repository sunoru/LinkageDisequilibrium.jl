"""
LD-Select Algorithm.

`haplotypes` is a nxm matrix (or n-element list of m-element array)
where the j-th element of the i-th row indicates the allele type of the j-th SNP
in the i-th haplotype.

Assuming that all the SNPs are biallelic, it does not matter whether 0's indicate MAF.

`MAF_threshold` is the minor-allele frequency threshold.

`r²_threshold` is the r² threshold.

Return:

`bins`: A dictionary from the indices of tagSNPs to the indices in that bin.
"""
function ld_select(
    haplotypes::Matrix{Bool};
    MAF_threshold = 0.1,
    r²_threshold = 0.8
)
    @assert r²_threshold < 1
    
    fs, haplotypes = prepare_select(haplotypes)
    n, m = size(haplotypes)
    snp_indices = (1:m)[fs .< (1 - MAF_threshold)]
    m2 = length(snp_indices)
    
    # R²[i,j] means r² between i-th and j-th SNPs > threshold.
    R² = r²(haplotypes[:, snp_indices])
    bR² = R² .> r²_threshold

    bins = Dict{Int, Vector{Int}}()
    unbinned = collect(1:m2)

    while length(unbinned) > 0
        tR² = bR²[unbinned, unbinned]
        p = vec(sum(tR², dims = 1))
        tagi = argmax(p)
        tag = snp_indices[unbinned[tagi]]
        binned_indices = vec(tR²[:, tagi])
        bins[tag] = unbinned[binned_indices]
        deleteat!(unbinned, binned_indices)
        p[tagi] ≡ 0 && break
    end
    for each in unbinned
        tag = snp_indices[each]
        bins[tag] = [tag]
    end
    bins
end

ld_select(
    haplotypes::AbstractVector;
    kwargs...
) = ld_select(make_binary(haplotypes); kwargs...)
