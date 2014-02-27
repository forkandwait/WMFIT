function [p D] = mkseries(T, iters=2000, shift=[], wiggle=0.1,  x0=[])
    ## mkseries(iters=2000, shift=[], wiggle=0.1,  x0=[])
    D = struct();
    
    ## initiate population output
    p = nan(rows(T), iters); 
    p(:,1) = ifelse(isempty(x0), rand(rows(T),1), x0);

    ## set up matrix for randomness
    foo = randn(rows(T), iters);

    ## set up matrix for mean-shift
    shift  = ifelse(isempty(shift), zeros(rows(T),1), shift);
    wiggle = ifelse(isempty(wiggle), 0, wiggle);


    ## fill population output 
    for ii = 1:(iters-1)
        p(:,ii+1) = T * p(:,ii) + foo(:,ii+1)  * wiggle + shift;
    end    
endfunction
