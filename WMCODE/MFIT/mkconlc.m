function out = mkconlc(speclist, matdims, x0=[], H=[], q=[], A_eqm=[], b_eqv=[], lb=[], ub=[], varargin)
    ## makes a bunch of constraint matrices and vectors based on a
    ## specification of the linear bases for parts of the matrix.
    ##
    ## Each spec is 
    ## 
    ##     spec.name = 'somename'                   ## just a string ID
    ##     spec.idx  = [1 2; 1 3; 1 4]              ## a N x 2 matrix of cell indices that are being constrainted in terms of T(n,m) (not T(:))
    ##     spec.vset = {
    ##          struct('lc', [5 3 2 1], 'vlb',  0.0, 'vub', 4.0),         ## a basis vector, variable lower bound, variable upper bound (multiplier)
    ##          struct('lc', [5 3 2 1], 'vlb',  0.0, 'vub', 4.0),        
    ##          struct('lc', [1 1 8 1], 'vlb', -1.0, 'vub', 1.0),
    ##                 }
    ##     spec.b    = [0 1 2];                       ## constant results for the equality matrix multiplication [I | p] * v = b. 
 
    ##  (TODO alternative way to set stuff: spec.vset(1).lc spec.vset.vlb spec.vset.vub)
    ## 
    ## "Speclist" is a 1 x M cell array of these specs
    ## 
    ## "matdims" is a 1 x 2 matrix = size(T) where T is the fitted
    ## transition matrix
    ## 
    ## "H" etc are the optional parameters to qp, which will be
    ## expanded to work with the additional slack variables
    ## 
    ## varargin collects whatever 

       
    ## format and error check

    ## determine size of everything that will be needed
    nslacks = 0;
    nvars   = 0;

    ## go through each spec
    for _spec = speclist(:)'
        spec = _spec{1};
        linmatrixthing = nan(0, columns(spec.vset{1}.lc));
        for  linthing = spec.vset(:)'
             nslacks += 1;      # count each vector component in a constraint set
             linmatrixthing(end+1,:) = linthing{1}.lc;
        endfor
        nvars   += rows(spec.idx);

        ## Check for stuff.  Difficult...
        if length(spec.b(:)) ~= rows(spec.idx)
           error("mkconlc: bad size for idx and b");
        endif
        if rank(linmatrixthing) < rows(linmatrixthing)
            warning("mkconlc: complete linear matrix of equality constraints is not full rank. Spec name: %s", spec.name);
        endif

    endfor

    ## allocate big matrix -- each cell in T gets a row, each slack gets a column
    nT = matdims(1)*matdims(2); # number of entries in transition matrix
    eqmat = zeros(nvars, nT + nslacks); 

    ## allocate constant vector and bounds 
    qvec  = nan(nvars,   1);
    lbvec = nan(nslacks, 1);
    ubvec = nan(nslacks, 1);


    ## take person-friendly input and make qp-friendly output.
    ## Pain in the ass!  One thing that helped was using separate counter 
    spcnt   = 0;               # which constraint set are we on
    conscnt = 0;               # which constraint item on (counted globally)
    bcnt    = 0;               # how many vector bases are we on 
    colstrt = 0;               # where we start in the equality matrix
    colstp  = 0;               # where we stop
    for _spec = reshape(speclist, 1, [])
        spec  = _spec{1};
        spcnt += 1; 

        lincm = [];
        for ii = 1:length(spec.vset)
            bcnt += 1;

            ## assemble LC components for easy assignment 
            lincm(:,ii) = (spec.vset{ii}.lc)(:);

            ## assign bounds
            lbvec(bcnt) = spec.vset{ii}.vlb;
            ubvec(bcnt) = spec.vset{ii}.vub; 

        endfor

        ## figure out where columns get placed
        colstrt = ifelse(spcnt == 1, nT+1, colstp+1);
        colstp  = colstrt + columns(lincm) - 1;

        if rows(spec.idx) ~= rows(lincm)
           error("mkconlc: bad spec dimensions. %i.", spcnt);
        endif 

        for ii = 1:rows(spec.idx);
            conscnt += 1;

            ## assign ones to appropriate places in eqmat for each T entry
            ti = sub2ind(matdims, spec.idx(ii,1), spec.idx(ii,2));
            eqmat(conscnt, ti) = 1;

            ## assign LC components down the right block 
            eqmat(conscnt, colstrt:colstp) = -1*lincm(ii, :); 

            ## assign constant term 
            qvec(conscnt) = spec.b(ii);

        endfor 
    endfor

    ## fix up some parameters for qp(), or their dimensions get
    ## screwed up. 
    out.x0     = vertcat(x0, rand(nslacks, 1) - 0.5);
    out.H      = blkdiag(H, zeros(nslacks));
    out.q      = vertcat(q, zeros(nslacks, 1));
    out.lbvec  = vertcat(lb, lbvec);
    out.ubvec  = vertcat(ub, ubvec);
    if ~isempty(A_eqm)
        A_eqm      = [A_eqm, zeros(rows(A_eqm), nslacks)];
        out.A_eqm  = vertcat(eqmat, A_eqm);
        out.b_eqv  = vertcat(qvec, b_eqv);
    else 
        out.A_eqm  = eqmat;
        out.b_eqv  = qvec;
    endif 


endfunction




%!shared spec1, spec2, spec3, spec4, spec5

%!     spec1      = struct();
%!     spec1.name = 'spec1';            
%!     spec1.idx  = [1 1; 1 2; 1 3; 1 4];    # which cells in the transition matrix
%!     spec1.vset = {                        
%!           struct('lc', [5 3 2 1], 'vlb',  0.0, 'vub', 4.0),       
%!           struct('lc', [1 1 1 1], 'vlb',  0.0, 'vub', 4.0),       
%!           struct('lc', [1 1 1 10], 'vlb', -1.0, 'vub', 1.0)
%!                 };
%!     spec1.b    = [0 1 2 0];    

%!     spec2      = struct();
%!     spec2.name = 'spec2';            
%!     spec2.idx  = [4 4];           
%!     spec2.vset = {struct('lc', 1, 'vlb', -1.0, 'vub', 1.0) };
%!     spec2.b    = 1;            

%!     spec3      = struct();
%!     spec3.name = 'spec3';            
%!     spec3.idx  = [ 1 2; 1 3];           
%!     spec3.vset = {struct('lc', [.5 .5], 'vlb', -1.0, 'vub', 1.0)};
%!     spec3.b    = [0 0];            


%!test
%!     foo = mkconlc({spec1}, [4 4]);

%!test
%!     foo = mkconlc({spec2}, [4 4]);

%!test
%!     foo = mkconlc({spec3}, [4 4]);

%!test     
%!     foo = mkconlc({spec1 spec2 spec3}, [4 4]);
