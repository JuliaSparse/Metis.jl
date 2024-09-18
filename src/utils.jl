# Only for square matrices. Metis 
function csr_to_csc(G)
    # We can get the transpose of G by assume its "CSR"
    Growval = G.rowval
    Growptr = G.colptr

    N = length(Growptr) - 1
    nnzG = length(Growval)

    GTcolptr = ones(idx_t,N + 1)
    GTrowval = zeros(idx_t,nnzG)

    # Count elements in each column:
    cnt = zeros(idx_t,N)
    for k in 1:nnzG
        col = Growval[k]
        cnt[col] += 1
    end
    # Cumulative sum to set the column pointer of matrix B
    for i in 2:N+1
        GTcolptr[i] = GTcolptr[i-1] + cnt[i-1]
    end

    for row in 1:N
        for j in (Growptr[row]:Growptr[row+1] - 1)
            col = Growval[j]
            dest = GTcolptr[col]

            GTrowval[dest] = row
            GTcolptr[col] += 1
        end
    end

    pop!(GTcolptr)
    GTcolptr = [idx_t(1); GTcolptr]

    return GTrowval, GTcolptr
end

function get_full_sparsity(G;uplo=:L)
    # Get the rowval and colptr from the triangle
    Growval = G.rowval
    Gcolptr = G.colptr
    # Get the rowval and colptr from the transpose of the triangle
    GTrowval, GTcolptr = csr_to_csc(G)

    # Allocating rowval and colptr for the full matrix (F)
    # Note that we might duplicate diagonal entries. This does not matter
    # as METIS ignores diagonal entries
    Frowval = zeros(idx_t, length(Growval) + length(GTrowval))
    Fcolptr = idx_t.(Gcolptr .- 1) + GTcolptr
    # Combining sparsity patterns:
    # Since we have to triangles we can combine by a "vcat"
    for col in 1:(length(Gcolptr) - 1)
        # For lower-triangular matrices the transpose comes first
        if uplo == :L
            ngt = GTcolptr[col + 1] - GTcolptr[col]
            for i = 1:ngt
                Frowval[Fcolptr[col] - 1 + i] = GTrowval[GTcolptr[col] + i - 1]
            end
            ng = Gcolptr[col + 1] - Gcolptr[col]
            for j = 1:ng
                Frowval[Fcolptr[col] - 1 + ngt + j] = Growval[Gcolptr[col] + j - 1]
            end
        # For upper-triangular matrices the transpose comes first
        else uplo == :U
            ng = Gcolptr[col + 1] - Gcolptr[col]
            for i = 1:ng
                Frowval[Fcolptr[col] - 1 + i] = Growval[Gcolptr[col] + i - 1]
            end
            ngt = GTcolptr[col + 1] - GTcolptr[col]
            for i = 1:ngt
                Frowval[Fcolptr[col] - 1 + ng + i] = GTrowval[GTcolptr[col] + i - 1]
            end
        end
    end

    return Frowval, Fcolptr

end
