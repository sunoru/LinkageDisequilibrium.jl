_ordered(a, b) = a ≤ b ? (a, b) : (b, a)
function _get_most_informative(d_edges, tagged_edges)
    s = 0
    si = 1
    for (i, edges) in enumerate(d_edges)
        t = sum([1 for e ∈ edges if e ∉ tagged_edges])
        if t > s
            s = t
            si = i
        end
    end
    si
end

"""
Selection with Informativeness Algorithm.

`haplotypes` is a nxm matrix (or n-element list of m-element array)
where the j-th element of the i-th row indicates the allele type of the j-th SNP
in the i-th haplotype.

Assuming that all the SNPs are biallelic, it does not matter whether 0's indicate MAF.

Note that currently this function implements a naive greedy algorithm,
which is expected to be time consuming.

`ξ` is the informativeness threshold.

Return:

`tagSNPs`: an array of SNP indices that has at least `ξ` informativeness.
"""
function informative_select(
    haplotypes::Matrix{Bool};
    MAF_threshold = 0.1,
    ξ = 0.8
)
    @assert ξ ≤ 1
    
    fs, haplotypes = prepare_select(haplotypes)
    n, m = size(haplotypes)
    snp_indices = (1:m)[fs .< (1 - MAF_threshold)]
    m2 = length(snp_indices)

    total = Set{Tuple{Int, Int}}()
    d_edges = Set{Tuple{Int, Int}}[]
    for index in snp_indices
        vs = haplotypes[:, index]
        v0s = (1:n)[vs .≡ false]
        v1s = (1:n)[vs .≡ true]
        edges = Set{Tuple{Int, Int}}()
        for i in v0s
            for j in v1s
                push!(edges, _ordered(i, j))
            end
        end
        push!(d_edges, edges)
        union!(total, edges)
    end
    len_d_edges = length.(d_edges)

    first = argmax(len_d_edges)

    tags = [snp_indices[first]]
    tagged_edges = d_edges[first]
    deleteat!(snp_indices, first)
    deleteat!(d_edges, first)

    target = ceil(Int, ξ * length(total))
    while length(tagged_edges) < target && length(tags) < m2
        si = _get_most_informative(d_edges, tagged_edges)
        union!(tagged_edges, d_edges[si])
        push!(tags, snp_indices[si])
        deleteat!(snp_indices, si)
        deleteat!(d_edges, si)
    end
    tags
end
 
informative_select(
    haplotypes::AbstractVector;
    kwargs...
) = informative_select(make_binary(haplotypes); kwargs...)
