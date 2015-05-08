module LDPCsim

using PyPlot

export
    ldpcMatrix,
    parityCheckMatrix,
    SymErr,
    sumProduct,
    sumProductLog,
    sumProductLogSimple,
    sumProductHard,
    ldpcEncode,
    ldpcSim


#########################################################################
function SymErr(x,y)
    # compare x and y to get # of symbol errors (bit errors for binary data)
    number = length(find(x.!=y));
    ratio= number/length(x);
    return (number,ratio)
end

#########################################################################
function ldpcMatrix(M, N, method, no4cycle, checksPerCol)
    # Create rate 1/2 LDPC adjacency matrix
    #
    # Inputs
    #  M        : num of rows
    #  N        : num of column
    #  method   : technique for distributing non-zero element
    #             0 evencol : For each column, place 1s uniformly at random
    #             1 evenboth: For each column and row, place 1s uniformly at random
    #  no4cyle  : if 1, then eliminate length-4 cycles, otherwise do nothing
    #  checksPerCol: num of ones per column
    #
    # Output
    #  H        : LDPC adjacency matrix                   

    # num of ones per row (N/M ratio must be 2)
    if N/M != 2
        @printf("Code rate must be 1/2 for now\n");
    end
    ## only used for evenboth
    checksPerRow = int((N/M)*checksPerCol);

    onesInCol = zeros(Int,M,N)
    if method==0   # evencol
        # Distribute 1's uniformly at random within column
        for nn = 1:N
            onesInCol[:, nn] = randperm(M);
        end

        # Create non zero elements (1s) index
        r = reshape(onesInCol[1:checksPerCol, :], N*checksPerCol, 1);
        tmp = repmat([1:N], checksPerCol, 1);
        c = reshape(tmp, N*checksPerCol, 1);

        # Create sparse matrix H
        H = full(sparse(r[:], c[:], ones(Int,N*checksPerCol), M, N,|));

    elseif method==1     # evenboth
        # Distribute 1s uniformly at random within column
        for nn = 1:N
            onesInCol[:, nn] = randperm(M)';
        end
        
        # Create non zero elements (1s) index
        r = reshape(onesInCol[1:checksPerCol, :], N*checksPerCol, 1);
        tmp = repmat([1:N], checksPerCol, 1);
        c = reshape(tmp, N*checksPerCol, 1);
        
        # Make the number of 1s between rows as uniform as possible     

        # Order row index
        ix = sortperm(r[:]);

        # Order column index based on row index
        cSort = zeros(Int,N*checksPerCol,1)
        for nn = 1:N*checksPerCol
            cSort[nn] = c[ix[nn]];
        end

        # Create new row index with uniform weight
# println("checksPerRow=$(checksPerRow)")
        tmp = repmat([1:M], checksPerRow, 1);
# println("tmp=$(tmp)")
# println("size(tmp)=$(size(tmp))")
# println("N=$(N)")
        r = reshape(tmp, N*checksPerCol, 1);

        # Create  H
        allOnes = ones(Int,N*checksPerCol,1)
        S = sparse(r[:], cSort[:], allOnes[:], M, N,|);
        H = full(S); 

    else
        error("unknown method $(method)")
    end

    # Check rows that have no 1 or only have one 1
    for mm = 1:M
        
        n = randperm(N);
        # Add two 1s if row has no 1
        if length(find(r == mm)) == 0
            H[mm, n[1]] = 1;
            H[mm, n[2]] = 1;
            # Add one 1 if row has only one 1   
        elseif length(find(r == mm)) == 1
            H[mm, n[1]] = 1;
        end

    end # for mm


    if no4cycle == 1  # if no4cycle, remove any length-4 cycle
        
        for mm = 1:M

            # Look for pair of row - column
            for jx = (mm + 1):M
                w = (&)(H[mm, :], H[jx, :]);
                c1 = find(w);
                lc = length(c1);
                if lc > 1
                    
                    # If found, flip one 1 to 0 in the row with less number of 1s
                    if length(find(H[mm, :])) < length(find(H[jx, :]))
                        # Repeat the process until only one column left 
                        for cc = 1:lc - 1
                            H[jx, c1[cc]] = 0;
                        end
                    else
                        for cc = 1:lc - 1
                            H[mm, c1[cc]] = 0;
                        end
                    end
                end
            end # for jx
        end # for mm
    end # if no4cycle==1

    return H
end

#########################################################################
function parityCheckMatrix(src, Hin, strategy)
    # Make parity check matrix from adjacency matrix H using LU
    # decomposition.
    #
    # Inputs
    #  src     : Binary source
    #  H       : adjacency matrix
    #  strategy: Strategy for finding the next non-zero diagonal elements
    #            0: First non-zero found by column search
    #            1: Minimum number of non-zeros in later columns
    #            2: Minimum product of:
    #               - num of non-zeros its column minus 1
    #               - num of non-zeros its row minus 1
    # Outputs
    #  c       : Check bits
    #  newH    : Rearrange H

    H = copy(Hin); # copy so as to not overwrite H

    # Get the matric dimension
    (M, N) = size(H);
    # Set a new matrix F for LU decomposition
    F = copy(H);
    # LU matrices
    L = zeros(M, N - M);
    U = zeros(M, N - M);

    # Re-order the M x (N - M) submatrix
    for mx = 1:M

        # strategy {0 = First; 1 = Mincol; 2 = Minprod}
        if strategy==0 # Create diagonal matrix using 'First' strategy
            
            # Find non-zero elements (1s) for the diagonal
            (r, c) = findn(F[:, mx:end]);
            
            # Find non-zero diagonal element candidates
            rowIndex = findn(r .== mx)[1];
            
            # Find the first non-zero column
            chosenCol = c[rowIndex[1]] + (mx - 1);

        elseif strategy==1 # Create diag matrix using 'Mincol' strategy
            
            # Find non-zero elements (1s) for the diagonal
            (r, c) = findn(F[:, mx:end]);
            colWeight = sum(F[:, mx:end], 1);
            
            # Find non-zero diagonal element candidates
            rowIndex = findn(r .== mx)[1];
            
            # Find the minimum column weight
            (x, ix) = findmin(colWeight[c[rowIndex]]);
            # Add offset to the chosen row index to match the dimension of the... 
            # original matrix F
            chosenCol = c[rowIndex[ix]] + (mx - 1);
            
        elseif strategy==2 # Create diag matrix using 'Minprod' strategy   
            
            # Find non-zero elements of F
            (r, c) = findn(F[:, mx:end]);
            colWeight = sum(F[:, mx:end], 1) - 1;
            rowWeight = sum(F[mx, :], 2) - 1;
            
            # Find non-zero diagonal element candidates
            rowIndex = findn(r .== mx)[1];
# if isempty(rowIndex)
#     println("M=$(M)")
#     println("mx=$(mx)")
#     println("r=$(r)")
# end
            # Find the minimum product
            (x, ix) = findmin(colWeight[c[rowIndex]]*rowWeight);
            # Add offset to the chosen row index to match the dimension of the
            # original matrix F
            chosenCol = c[rowIndex[ix]] + (mx - 1);
            
        else
            @printf("Please select columns re-ordering strategy!\n");
            
        end # switch

        # Re-ordering columns of both H and F
        tmp1 = F[:, mx];
        tmp2 = H[:, mx];
        F[:, mx] = F[:, chosenCol];
        H[:, mx] = H[:, chosenCol];
        F[:, chosenCol] = tmp1;
        H[:, chosenCol] = tmp2;
        
        # Fill the LU matrices column by column
# println("mx=$(mx)")
# println("size(mx)=$(size(mx))")
# println("L=$(L)")
# println("size(L)=$(size(L))")
# println("F=$(F)")
# println("size(F)=$(size(F))")
# println("U=$(U)")
# println("size(U)=$(size(U))")
        L[mx:end, mx] = F[mx:end, mx];
        U[1:mx, mx] = F[1:mx, mx];
        
        # There will be no rows operation at the last row
        if mx < M           
            
            # Find the later rows with non-zero elements in column mx
            r2 = findn(F[(mx + 1):end, mx])[1]; 
            # Add current row to the later rows which have a 1 in column mx
            F[(mx + r2), :] = mod(F[(mx + r2), :] + repmat(F[mx, :], length(r2), 1), 2);
            
        end # if
        
    end # for mx
# println("H=$(H)")
# println("L=$(L)")
# println("U=$(U)")
# println("L*U=$(L*U)")
    # Find src
    z = mod(H[:, (N - M) + 1:end]*src, 2);

    # Parity check vector found by solving sparse LU
    c = mod(U\(L\z), 2); 

    # Return the rearrange H 
    newH = H;

    return (c, newH)
end


#########################################################################
function sumProduct(zdata, H, N0, Niter)
    # Probability-domain sum product algorithm LDPC decoder
    #
    # Inputs
    #   zdata     : Received signal vector (column vector)
    #   H         : LDPC matrix
    #   N0        : Noise variance
    #   Niter     : num of iterations
    #
    # Outputs
    #   uhat      : Decoded vector (0/1) 
    (M,N) = size(H);

    # priors
    P1 = ones(size(zdata))./(1 + exp(-2*zdata./(N0/2)));
    P0 = 1 - P1;

    # Init
    K0 = zeros(M, N);
    K1 = zeros(M, N);
    R0 = zeros(M, N);
    R1 = zeros(M, N);
    Q0 = H.*repmat(P0', M, 1);
    Q1 = H.*repmat(P1', M, 1);
    uhat = zeros(N,1);

    for n = 1:Niter
        
        # index over columns
        for mx = 1:M
            
            # Find non-zeros in the column
            c1 = find(H[mx, :]);

            for k = 1:length(c1)
                
                # Get column products of dR\c1(l)
                dR = 1.0;
                for l = 1:length(c1)
                    if l!= k
                        dR = dR*(Q0[mx, c1[l]] - Q1[mx, c1[l]]);
                    end
                end # for l
                
                R0[mx, c1[k]] = (1 + dR)/2;
                R1[mx, c1[k]] = (1 - dR)/2;
                
            end # for k
            
        end # for mx
        
        # index over rows
        for nx = 1:N
            
            # Find non-zeros in the row
            r1 = find(H[:, nx]);
            
            for k = 1:length(r1)
                
                # Get row products of prodOfrij\ri(l)
                prodOfrij0 = 1;
                prodOfrij1 = 1;   
                for l = 1:length(r1)
                    if l!= k
                        prodOfrij0 = prodOfrij0*R0[r1[l], nx];
                        prodOfrij1 = prodOfrij1*R1[r1[l], nx];
                    end
                end # for l
                
                # Update constants
                K0[r1[k], nx] = P0[nx]*prodOfrij0;
                K1[r1[k], nx] = P1[nx]*prodOfrij1;
                
                # Update Q0 and Q1
                Q0[r1[k], nx] = K0[r1[k], nx]./(K0[r1[k], nx] + K1[r1[k], nx]);
                Q1[r1[k], nx] = K1[r1[k], nx]./(K0[r1[k], nx] + K1[r1[k], nx]);
                
            end # for k
            
            # Update constants
            Ki0 = P0[nx]*prod(R0[r1, nx]);
            Ki1 = P1[nx]*prod(R1[r1, nx]);
            
            # Get Qj
            Qi0 = Ki0/(Ki0 + Ki1);
            Qi1 = Ki1/(Ki0 + Ki1);
            
            # Decode Qj        
            if Qi1 > Qi0
                uhat[nx] = 1;
            else
                uhat[nx] = 0;
            end
            
        end # for nx
        
    end # for n

    return uhat
end


#########################################################################
function sumProductLog(zdata, H, N0, Niter)
    # Log-domain sum product algorithm LDPC decoder
    #
    # Inputs
    #  zdata     : Received signal vector (column vector)
    #  H         : LDPC matrix
    #  N0        : Noise variance
    #  Niter     : num of iterations
    #
    # Outputs
    #  uhat      : Decoded vector (0/1) 

    (M,N) = size(H);

    # prior log-likelihood. Minus sign is used for 0/1 to -1/1 mapping
    Lci = (-4*zdata./N0)';

    # Init
    LR = zeros(M, N);
    Pibetaij = zeros(M, N);
    uhat = zeros(N,1);

    # Asscociate the L(ci) matrix with non-zero elements of H
    LQ = H.*repmat(Lci, M, 1);
    
    # Get non-zero elements
    (r, c) = findn(H);

    for n = 1:Niter
        
        # Get the sign and magnitude of L(Q)   
        alphaij = sign(LQ);   
        betaij = abs(LQ);

        for l = 1:length(r)
            Pibetaij[r[l], c[l]] = log((exp(betaij[r[l], c[l]]) + 1)/
                                       (exp(betaij[r[l], c[l]]) - 1));
        end
        
        # index over columns
        for i = 1:M
            
            # Find non-zeros in the column
            c1 = find(H[i, :]);

            # Get the summation of Pi(betaij))        
            for k = 1:length(c1)

                sumOfPibetaij = 0;
                prodOfalphaij = 1;

                # Summation of Pi(betaij)\c1[k]
                Pic1 = Pibetaij[i, c1]
                sumOfPibetaij = sum(Pic1) - Pic1[1,k];
#                sumOfPibetaij = sum(Pibetaij[i, c1]) - Pibetaij[i, c1[k]];
                
                # Avoid division by zero/very small number, get Pi(sum(Pi[betaij]))
                if sumOfPibetaij < 1e-20
                    sumOfPibetaij = 1e-10;
                end         
                PiSumOfPibetaij = log((exp(sumOfPibetaij)+1)/(exp(sumOfPibetaij) - 1));
                
                # Multiplication of alphaij\c1[k] (use '*' since alphaij are -1/1s)
                prodOfalphaij = prod(alphaij[i, c1])*alphaij[i, c1[k]];
                
                # Update L[R]
                LR[i, c1[k]] = prodOfalphaij*PiSumOfPibetaij;
                
            end # for k
            
        end # for i

        # index over rows
        for j = 1:N

            # Find non-zero in the row
            r1 = find(H[:, j]);
            
            for k = 1:length(r1)        
                
                # Update L[Q] by summation of L[rij]\r1[k]
                LQ[r1[k], j] = Lci[j] + sum(LR[r1, j]) - LR[r1[k], j];
                
            end # for k
            
            # Get L[Qi]
            LQi = Lci[j] + sum(LR[r1, j]);
            
            # Decode L(Qi)
            if LQi < 0
                uhat[j] = 1;
            else
                uhat[j] = 0;
            end
            
        end # for j
        
    end # for n

    return uhat
end

#########################################################################
function sumProductLogSimple(zdata, H, Niter)
    # Simplified log-domain sum product algorithm LDPC decoder
    #
    # Inputs
    #  zdata     : Received signal vector (column vector)
    #  H         : LDPC matrix
    #  Niter     : num of iterations
    #
    # Outputs
    #  uhat      : Decoded vector (0/1)

    (M,N) = size(H);

    # Prior log-likelihood (simplified). Minus sign is used for 0/1 to -1/1 mapping
    Lci = -zdata';

    # Init
    LR = zeros(M, N);
    Pibetaij = zeros(M, N);
    uhat = zeros(N,1);

    # Asscociate the L(ci) matrix with non-zero elements of H
    LQ = H.*repmat(Lci, M, 1);

    for n = 1:Niter
        
        
        # Get the sign and magnitude of L(Q)   
        alphaij = sign(LQ);   
        betaij = abs(LQ);

        # index over columns
        for i = 1:M
            
            # Find non-zeros in the column
            c1 = find(H[i, :]);
            
            # Get the minimum of betaij
            for k = 1:length(c1)
                
                # Minimum of betaij\c1[k]
                minOfbetaij = realmax();
                for l = 1:length(c1)
                    if l != k 
                        if betaij[i, c1[l]] < minOfbetaij
                            minOfbetaij = betaij[i, c1[l]];
                        end
                    end          
                end # for l
                
                # Multiplication alphaij\c1(k) (use '*' since alphaij are -1/1s)
                prodOfalphaij = prod(alphaij[i, c1])*alphaij[i, c1[k]];
                
                # Update L[R]
                LR[i, c1[k]] = prodOfalphaij*minOfbetaij;
                
            end # for k
            
        end # for i

        # index over rows
        for j = 1:N

            # Find non-zero in the row
            r1 = find(H[:, j]);
            
            for k = 1:length(r1)        
                
                # Update L[Q] by summation of L[rij]\r1[k]
                LQ[r1[k], j] = Lci[j] + sum(LR[r1, j]) - LR[r1[k], j];
                
            end # for k
            
            # Get L[Qi]
            LQi = Lci[j] + sum(LR[r1, j]);
            
            # Decode L(Qi)
            if LQi < 0
                uhat[j] = 1;
            else
                uhat[j] = 0;
            end
            
        end # for j
        
    end # for n

    return uhat
end

#########################################################################
function sumProductHard(zdata, H, Niter)
    # Hard-decision/bit flipping sum product algorithm LDPC decoder
    #
    # Inputs
    #  zdata     : Received signal vector (column vector)
    #  H         : LDPC matrix
    #  Niter     : num of iterations
    #
    # Outputs
    #  uhat      : Decoded vector (0/1) 

    (M,N) = size(H);

    # Prior hard-decision
    ci = 0.5*(sign(zdata') + 1);

    # Init
    R = zeros(M, N);
    uhat = zeros(N,1);

    # Asscociate the ci matrix with non-zero elements of H
    Q = H.*repmat(ci, M, 1);
    
    for n = 1:Niter
        
        
        # index over columns
        for i = 1:M
            
            # Find non-zeros in the column
            c1 = find(H[i, :]);
            
            # Get the summation of Q\c1(k)        
            for k = 1:length(c1)

                R[i, c1[k]] = mod(sum(Q[i, c1]) + Q[i, c1[k]], 2);
                
            end # for k
            
        end # for i
        
        # index over rows
        for j = 1:N

            # Find non-zero in the row
            r1 = find(H[:, j]);

            # Number of 1s in a row
            numOfOnes = length(find(R[r1, j]));

            for k = 1:length(r1) 

                # Update Q, set '1' for majority of 1s else '0', excluding r1(k)
                if numOfOnes + ci[j] >= length(r1) - numOfOnes + R[r1[k], j]
                    Q[r1[k], j] = 1;
                else
                    Q[r1[k], j] = 0;
                end
                
            end # for k
            
            # Bit decoding
            if numOfOnes + ci[j] >= length(r1) - numOfOnes
                uhat[j] = 1;
            else
                uhat[j] = 0;
            end
            
        end # for j
        
    end # for n
    return uhat
end

#########################################################################
# @doc doc"""ldpcSim(Ns,rho, M,Niter,rate,Mc,mSparse,strategy)
# Simulate transmission of LDPC-coded signals over an AWGN channel
#
# Examples show SER as we deviate away from LTE PUSCH defaults
#  # Performance as a function of SNR
#  ldpcSim(400,[-5.0:5:20],[10])
#  # Performance as a function of M
#  ldpcSim(1000,[3.5],[4,8,16,32,64,128,256])
#  # Performance as a function of SNR for two different M values
#  ldpcSim(400,[-5.0:5:20],[10,100])
# """ ->
function ldpcSim{Ti<:Integer,Tf<:FloatingPoint}(Ns::Ti,rho::Array{Tf,1},
                M::Array{Ti,1};Niter::Ti=5,rate::Tf=.5,Mc::Ti=2,mSparse::Ti=0,
                strategy::Ti=2,Nmat=1)

    println("   Ns   rho     M   Niter rate    Mc mSparse strategy")

    out = cell(0)
    if length(rho)>1 && length(M)==1 && Nmat==1
        for s = 1:length(rho)
            push!(out,ldpcEach(Ns,Mc,M[1],rate,mSparse,strategy,Niter,rho[s]));
        end
        xval = rho;
        xlab = "rho (Eb/N0 in dB)";
        tstr = @sprintf("Mc=%d,Ns=%d,M=%d,rate=%3.2f,mSparse=%d,strategy=%d,Niter=%d",
                        Mc,Ns,M[1],rate,mSparse,strategy,Niter);
    elseif length(rho)>1 && length(M)>1 && Nmat==1
        for mx = 1:length(M)
            for rx = 1:length(rho)
                outi=ldpcEach(Ns,Mc,M[mx],rate,mSparse,strategy,Niter,rho[rx]);
                for a=1:length(outi[2])
                    outi[2][a] = outi[2][a] * ", M=$(M[mx])"
                end
                if mx==1
                    push!(out,outi)
                else
                    t1 = [out[rx][1]; outi[1]]
                    t2 = [out[rx][2]; outi[2]]
                    out[rx] = (t1,t2)
                end
            end
        end
        xval = rho;
        xlab = "rho (Eb/N0 in dB)";
        tstr = @sprintf("Mc=%d,Ns=%d,rate=%3.2f,mSparse=%d,strategy=%d,Niter=%d",
                        Mc,Ns,rate,mSparse,strategy,Niter);
    elseif length(M)>1
        for s = 1:length(M)
            push!(out,ldpcEach(Ns,Mc,M[s],rate,mSparse,strategy,Niter,rho[1]));
        end
        xval = M;
        xlab = "M: message size";
        tstr = @sprintf("Mc=%d,Ns=%d,rate=%3.2f,mSpar=%d,strat=%d,Niter=%d,Eb/N0=%4.1f",
                        Mc,Ns,rate,mSparse,strategy,Niter,rho[1]);
    elseif length(rate)>1
        for s = 1:length(rate)
            push!(out,ldpcEach(Ns,Mc,M[1],rate[s],mSparse,strategy,Niter,rho[1]));
        end
        xval = rate;
        xlab = "rate: M/N";
        tstr = @sprintf("Mc=%d,Ns=%d,M=%d,mSparse=%d,strategy=%d,Niter=%d,Eb/N0=%4.1f",
                        Mc,Ns,M[1],mSparse,strategy,Niter,rho[1]);
    elseif length(Niter)>1
        for s = 1:length(Niter)
            push!(out,ldpcEach(Ns,Mc,M[1],rate,mSparse,strategy,Niter[s],rho[1]));
        end
        xval = Niter;
        xlab = "Niter";
        tstr = @sprintf("Mc=%d,Ns=%d,M=%d,mSparse=%d,strategy=%d,rate=%3.2f,Eb/N0=%4.1f",
                        Mc,Ns,M[1],mSparse,strategy,rate,rho[1]);
    elseif length(Mc)>1
        for s = 1:length(Mc)
            push!(out,ldpcEach(Ns,Mc[s],M[1],rate,mSparse,strategy,Niter,rho[1]));
        end
        xval = log2(Mc);
        xlab = "# constellation bits";
        tstr = "";
    elseif Nmat>1
        for nm = 1:Nmat
            srand(nm)
            for s = 1:length(rho)
                outi = ldpcEach(Ns,Mc,M[1],rate,mSparse,strategy,Niter,rho[s]);
                for a=1:length(outi[2])
                    outi[2][a] = outi[2][a] * " sd$(nm)"
                end
                if nm==1
                    push!(out,outi)
                else
                    for a=1:length(outi[2])
                        t1 = [out[s][1]; outi[1]]
                        t2 = [out[s][2]; outi[2]]
                        out[s] = (t1,t2)
                    end
                end
            end
        end
        xval = rho;
        xlab = "rho (Eb/N0 in dB)";
        tstr = @sprintf("Mc=%d,Ns=%d,M=%d,rate=%3.2f,mSparse=%d,strategy=%d,Niter=%d",
                        Mc,Ns,M[1],rate,mSparse,strategy,Niter);
    else
        @printf("At least one argument must be a vector")
    end

    pColor = {"r>-", "bo--","kx-.","gd-", "c^--","m*-.",
              "rs--","gp-.","bv-", "kh--","c+-.","m.-",};
    pIdx = 1;

    # Find number of samples
    Ns = size(out,1);

    #clf()
    ser = zeros(Ns);
    for a=1:length(out[1][2])
        for k=1:Ns
            ser[k] = out[k][1][a];
        end
        semilogy(xval,ser,pColor[pIdx],label=out[1][2][a]);
        hold(true);
        if pIdx == length(pColor)
            pIdx = 1;
        else
            pIdx = pIdx + 1;
        end
    end

    hold(false);
    xlabel(xlab);
    ylabel("BER");
    legend(loc=0);
    grid(which="both",axis="y")
    grid(which="major",axis="x")
    title(tstr)

    return
end



#######################################################################
function ldpcEach(Ns,Mc,M,rate,mSparse,strategy,Niter,rdb)

    rho     = 10^(rdb/10);

    if rate != 2.0
        @sprintf("Only rate 1/2 tested so far.")
    end

    N = integer(ceil(M/rate));

    # mSparse is method for sparse matrix creation (0 = evencol; 1 = evenboth)
    # Niter is number of iterations
    # strategy is matrix reorder strategy (0 = 1st; 1 = Mincol; 2 = Minprod)
    no4cycle = 1; # if 1, remove length-4 cycles
    checksPerCol = 3; # Number of 1s per column for LDPC matrix

    # # Make the LDPC matrix
    # H = ldpcMatrix(M, N, mSparse, no4cycle, checksPerCol);

    Nalgs = 3;
    ber = zeros(Nalgs);
    times = zeros(Ns);
    algNames = cell(Nalgs)
    u = zeros(N)
    newH = zeros(M,N);
    #Profile.clear()

    if isinteractive(); @printf("     "); end

    # Make random data (0/1)
    dSource = round(rand(M, Ns));
    for jx = 1:Ns
        tic();
        if (Ns<200 || mod(jx,100)==0) && isinteractive()
            @printf("\b\b\b\b\b%5d", jx);
        end

        # Ugly temporary trick to handle parity check errors; 
        flag=0
        while flag>-1
            try
                (u,newH) = ldpcEncode(M, N, mSparse, no4cycle,
                             checksPerCol,dSource[:, jx],strategy);
                flag=-1
            catch
                flag+=1
            end
        end
        # # Make the LDPC matrix
        # H = ldpcMatrix(M, N, mSparse, no4cycle, checksPerCol);
        # while rank(H)<M
        #     H = ldpcMatrix(M, N, mSparse, no4cycle, checksPerCol);
        # end
        # # Encoding message
        # (c, newH) = parityCheckMatrix(dSource[:, jx], H, strategy);
        # u = [c; dSource[:, jx]];

        # BPSK modulation
        bpskMod = 2*u - 1;
        # white gaussian noise
        N0 = 1/rho;
        tx = bpskMod + sqrt(N0)*randn(size(bpskMod));

        ax= 0;

        ax = ax+1;
        algNames[ax] = "uncoded";
        uhat = (sign(tx)+1)/2
        (num, rat) = SymErr(uhat[:], u); #  includes parity bits
        ber[ax] += rat;

        ax = ax+1;
        algNames[ax] = "sumProductLog";
        #algNames[ax] = "spl";
        uhat = sumProductLog(tx, newH, N0, Niter);
        (num, rat) = SymErr(uhat[:], u); #  includes parity bits
        ber[ax] += rat;

        ax = ax+1;
        algNames[ax] = "sumProductSimpler";
        #algNames[ax] = "sps";
        uhat = sumProductLogSimple(tx, newH, Niter);
        (num, rat) = SymErr(uhat[:], u);
        ber[ax] += rat;

        # ax = ax+1;
        # algNames[ax] = "sumProduct";
        # uhat = sumProduct(tx, newH, N0, Niter);
        # (num, rat) = SymErr(uhat[:], u);
        # ber[ax] += rat;

        # Better than uncoded only @ high SNR; after the threshold/waterfall
        # ax = ax+1;
        # algNames[ax] = "sumProductHard";
        # uhat = sumProductHard(tx, newH, Niter);
        # (num, rat) = SymErr(uhat[:], u);
        # ber[ax] += rat;

        times[jx] = toq();
        if sum(times)>Ns/40
            if isinteractive(); @printf("\b\b\b\b\b"); end
            Ns = jx
            break
        end
    end # for jx

    ber = ber/Ns;

    # using ProfileView
    # ProfileView.view()

    @printf("%5d %5.1f %5d %5d   %3.2f %5d %5d  %5d  ",
            Ns, rdb,  M, Niter, rate, Mc, mSparse, strategy)

    for ax = 1:Nalgs
        if ax<6
            @printf("    %7.4f",ber[ax])
        end 
    end
    @printf("\n")

    return ber, algNames
end

#######################################################################
function ldpcEncode(M, N, mSparse, no4cycle, checksPerCol,src,strategy)
    # Make the LDPC matrix
    H = ldpcMatrix(M, N, mSparse, no4cycle, checksPerCol);
    while rank(H)<M  # Hack to handle low-rank H.  
        H = ldpcMatrix(M, N, mSparse, no4cycle, checksPerCol);
    end
    # Encoding message
    (c, newH) = parityCheckMatrix(src, H, strategy);
    u = [c; src];
    return (u,newH)
end


#########################################################################
end # LDPCsim module
