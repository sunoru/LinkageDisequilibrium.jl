using Distributions

function fisher_exact_test(contingency_table)
    @assert size(contingency_table) ≡ (2, 2)
    w, y, x, z = contingency_table
    n = w + y + x + z
    p = binomial(w + x, w) * binomial(y + z, y) / binomial(n, w + y)
end

function χ²_test(contingency_table)
    @assert size(contingency_table) ≡ (2, 2)
    w, y, x, z = contingency_table
    n = w + y + x + z
    fA = (w + x) / n
    fB = (w + y) / n

    expected = [
        fA * fB        fA * (1 - fB)
        (1 - fA) * fB  (1 - fA) * (1 - fB)
    ] * n
    χ² = sum(contingency_table .^ 2 ./ expected) - n
    ccdf(Chisq(1), χ²)
end
