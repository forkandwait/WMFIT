 function [T V DEBUGV] = fitm(D0, D1,
                              eq_spec=[],
                              lc_spec=[],
                              lb_spec=[], 
                              ub_spec=[],
                              ineq_spec=[],
                              varargin)
     ## fit markov chains in a general way.
     ##
     ## Cell equality, cell inequality, cell linear combination
     ## equality. 
          
     ## check params (incomplete)
     if !all(size(D0) == size(D1))
         error("fitm: input and output data matrices (D0, D1) must have the same dimension (TODO: measurement mat)");
     end
     
    ## handle optional args
    if isempty(varargin)
        OPTS = struct();
    else
        OPTS = struct(varargin{:});
    endif

    ## median normalize
    if isfield(OPTS, 'median') && OPTS.median == 0 
        warning ("fitm: NOT using median normalization");
        d0 = D0;
        d1 = D1; 
    else 
        dmed = median([D0; D1](:));
        d0 = D0 ./ dmed;
        d1 = D1 ./ dmed; 
    endif

    ## use optional sqp approaches
    if isfield(OPTS, 'algo') && strmatch('sqp', OPTS.algo)
        1;
    else
        OPTS.algo = 'qp';
    endif

    ## set bicriterion regularization proportion -- defaults to all fit
    if isfield(OPTS, 'alpha')
        1;
    else
        OPTS.alpha = [];
    endif
    lccon.alpha = OPTS.alpha;

    ## Per-cell equality constraint thing, also gives sparsity.
    ## Remember to pass to mkconlc ..
    if ~isempty(eq_spec)

        ## make matrices out of spec
        eqcon = mkconeq(eq_spec);
        A_eqm = eqcon.m;
        b_eqv = eqcon.v;

        ## Determine which cells can be ignored due to sparsity. Note that sparse_r is an intermediate result
        sparse_r = logical(sum(abs(sign(A_eqm)),2) & (b_eqv == 0));
        sparse_c = logical(sum(A_eqm(sparse_r,:))'); 

        ## MODIFY A_eqm and b_eqv from sparsity (duh... tricky ...)
        A_eqm = A_eqm(!sparse_r,:);
        b_eqv = b_eqv(!sparse_r,:);
    else
        A_eqm = [];
        b_eqv = [];
        sparse_r = logical([]);
        sparse_c = logical(zeros(rows(d0),1));
    end

    ## handle initial guess for x0
    if isfield(OPTS, 'startx0') 
        x0 = OPTS.startx0(:);
    else
        x0 = zeros(rows(d0(!sparse_c,:))^2, 1); 
    endif 

    ## set up basic variables
    kI   = speye(rows(d0));
    M    = kron(d0', kI);  
    M_sp = M(:, !sparse_c);
    z    = sparse(d1(:));

    ## set up basic qp inputs with sparsity
    H = (M_sp'*M_sp);                 # note: qp() multiplies H by 0.5 so we don't here
    q = (-(z'*M_sp)');

    ## .. lower, upper bounds, forcing sparsity ..
    if isempty(lb_spec)
        lb = [];
    else 
        lb = lb_spec(:)(!sparse_c);
    endif
    if isempty(ub_spec)
        ub = [];
    else 
        ub = ub_spec(:)(!sparse_c);
    endif

    ## .. linear combination constraint thing.
    if ~isempty(lc_spec) && ~strcmpi(class(lc_spec), "cell")
        error("fitm: lc_spec must be a cell array with one or more linear combo constraints");
    end
    lccon = mkconlc(lc_spec, [rows(D0) rows(D0)], x0, H, q, A_eqm, b_eqv, lb, ub);
    lccon.q = sparse(lccon.q);
    lcnum = rows(lccon.q) - rows(q);
keyboard
    ## set up inequality matrices from spec.  SEEMS SHOULD BE APPENDED....
    if ~isempty(ineq_spec) && strmatch(OPTS.algo, 'qp')
       ineqcon = mkconineq(ineq_spec);
       lccon.A_in = ineqcon.m;
       lccon.A_lb = ineqcon.v1;
       lccon.A_ub = ineqcon.v2;
    elseif ~isempty(ineq_spec) && strmatch(OPTS.algo, 'sqp')
       ineqcon = mkconineq(ineq_spec, 1);
       lccon.A_in = ineqcon.m;
       lccon.A_lb = ineqcon.v;
    else
       lccon.A_in = [];            
       lccon.A_lb = [];            
       lccon.A_ub = [];            
    endif
        
keyboard

    ## do the qp, use return values
    [x D] = wmqp(OPTS.algo, lccon);
    if strmatch(OPTS.algo, 'qp')
        switch  D.info.info 
            case 0
                 ;
            case 1
                warning("fitm: qp not convex, local solution found, solveiter=%i", D.info.solveiter);
            case 2
                error(  "fitm: qp not convex and unbounded");
            case 3 
                warning("fitm: max iterations reached, solveiter=%i", D.info.solveiter);
            case 6
                error(  "fitm: infeasible qp");
            otherwise
                error(  "fitm: weird qp code: code=%i", D.info.info);
        endswitch
    elseif strmatch(OPTS.algo, 'sqp')
        switch  D.info 
            case 101
                 ;
            case 102
                error("fitm: sqp BFGS update failed");
            case 103
                warning("fitm: sqp max iterations reached");
            otherwise
                error("fitm: weird sqp code: code=%i", D.info);
        endswitch
    endif
           
    
    ## return stuff, including decompose and analyze obj and lambdas
    ## XXX DEAL WITH SPARSITY IN REBUILDING ... NOT YET...

    ## transition matrix, not including "special" variables
    _tdim = rows(d1)*rows(d0);
    T = reshape(x(1:_tdim), rows(d1), rows(d0));

    ## norm of fit calculated from x0

    ## special variables 
    if _tdim < length(x(:))
        V = x(:)(_tdim+1:end);
    else 
        V=[];
    endif

    ## diagnostics 
    _ssqD1     = sumsq(D1(:));
    _term1     = lccon.q' * x;
    _term2     = 0.5 * x'*lccon.H*x;
    _obj       = _term1 + _term2;
    _ssqD1mod  = sumsq((T*D0)(:));

    ## dump all variables if want
    if nargout >= 3
       DEBUGV = allv();
    endif
           
endfunction



## function [T DEBUGV] = fitm(D0, D1, eq_spec=[], lc_spec=[], ineq_spec=[], lb_spec=[], ub_spec=[], varargin) 


%!shared eq_spec1, lc_spec1, lb_spec1, ub_spec1, lb_spec2, ub_spec2, p1, T1, p2, T2
%!  eq_spec1 = {0, {}, {}, 0; {}, 0, 0, 0; 0, {}, 0, 0; 0, 0, {}, 0};

%!  lc_spec1 = struct();
%!  lc_spec1.name = 'fert';
%!  lc_spec1.idx  = [1 2; 1 3]; 
%!  lc_spec1.vset = { struct('lc', [.5 .5], 'vlb', .5, 'vub', 4) };
%!  lc_spec1.b    = [0 0];

%!  lb_spec1 = zeros(4);
%!  ub_spec1 = ones(4);

%!  lb_spec2 = repmat(-1, 4, 4);
%!  ub_spec2 = repmat(1, 4, 4);
  
%!test
%!  p = nan(4,6);
%!  p(:,1) = 1;
%!  T = [0 0.7 0.7 0; .95 0 0 0; 0 .95 0 0; 0 0 .5 0];
%!  for ii = 1:15
%!      p(:,ii+1) = T * p(:,ii);
%!  end
%!  [m D] = fitm(p(:,1:end-1), p(:,2:end), [], {lc_spec1}, [], lb_spec1, ub_spec1);



    
##%!test
##%!  d = rand(4,8);
##%!  [m D] = fitm(d, d);
##%!  assert(m*d, d, 0.000001);
## 
##%!test 
##%!  p = nan(4,6);
##%!  p(:,1) = 1;
##%!  T = [0 0 0 1; 1 0 0 0; 0 1 0 0; 0 0 1 0];
##%!  for ii = 1:20
##%!      p(:,ii+1) = T * p(:,ii) + (rand(4,1)-0.5)*0.1;
##%!  end
##%!  [m D] = fitm(p(:,1:end-1), p(:,2:end), eq_spec1, [], [], [], [], 'treg', 0.0, 'median', 0);
    
 
##%!test
##%!  p = nan(4,6);
##%!  p(:,1) = 1;
##%!  T = [0 0.7 0.7 0; .95 0 0 0; 0 .95 0 0; 0 0 .5 0];
##%!  for ii = 1:15
##%!      p(:,ii+1) = T * p(:,ii);
##%!  end
##%!  [m D] = fitm(p(:,1:end-1), p(:,2:end), [], [], [], lb_spec1, ub_spec1);
## 
## 
##%!test  
##%!  p = nan(4,6);
##%!  p(:,1) = .25;
##%!  T = [0 0.7 0.7 0; .95 0 0 0; 0 .95 0 0; 0 0 .5 0];
##%!  for ii = 1:15
##%!      p(:,ii+1) = T * p(:,ii);
##%!  end
##%!  [m D] = fitm(p(:,1:end-1), p(:,2:end), [], {lc_spec1}, [], lb_spec1, ub_spec1);
##%!  assert(m, T, .000001)
