#addpath (genpath('e:\WoodsMethod'));

## experiment with sqp regularized stuff
if 1

    ## basic data
    T  = [0.0000 0.0000 0.0143; 
          1.0000 0.4611 0.0000; 
          0.0000 0.5389 0.9664];
    rand("seed", 1);
    p  = mkseries(T, 5, [], [], [5.32 24.84 115.5]');

    ########################
    ## with no bugs, my code works sort of like his...
    eq_spec   = [nan 0 nan; nan nan 0; 0 nan nan];
    lb_spec   = zeros(3,3);
    
    ineq_spec = {
                 {[1 0 0; 1 0 0; 0 0 0], -10000, 1}, 
                 {[0 0 0; 0 1 0; 0 1 0], -10000, 1}, 
                 {[0 0 0; 0 0 0; 0 0 1], -10000, 1}
                 }; 

    ##                                                eq  lc       lb  ub       ineq
    ##[T_fit D]   = fitm(p(:,1:end-1), p(:,2:end), eq_spec, [], lb_spec, [], ineq_spec, 'algo', 'sqp1');

    [T1_fit D1] = fitm(p(:,1:end-1), p(:,2:end), eq_spec, [], lb_spec, [], ineq_spec, 'algo', 'sqp2', 'alpha', 1);
endif


## randomness
if 0
    T0        = diag([1 1 1], -1);
    T0(1,2:3) = 0.5;
    
    T1 = [ 
           0.00000   0.58375   0.58375   0.00000
           0.90000   0.00000   0.00000   0.00000
           0.00000   0.90000   0.00000   0.00000
           0.00000   0.00000   0.90000   0.00000
    ];
    T2 = [ 
           0.00000   0.63450   0.63450   0.00000
           0.90000   0.00000   0.00000   0.00000
           0.00000   0.75000   0.00000   0.00000
           0.00000   0.00000   0.50000   0.00000
    ];
    
    #eq_spec1 = {0, {}, {}, 0; {},  0, 0, 0; 0, {},  0,  0; 0, 0, {},  0};
    eq_spec1 = zeros(4);
    eq_spec1(1, 2:3) = nan;
    eq_spec1(2,1) = eq_spec1(3,2) = eq_spec1(4,3) = nan;
    
    lc_spec1 = struct();
    lc_spec1.name = 'fert';
    lc_spec1.idx  = [1 2; 1 3]; 
    lc_spec1.vset = { struct('lc', [.5 .5], 'vlb', .5, 'vub', 4) };
    lc_spec1.b    = [0 0];
    
    lb_spec1 = zeros(4);
    ub_spec1 = ones(4);
    
    lb_spec2 = repmat(-1, 4, 4);
    ub_spec2 = repmat(1, 4, 4);

    p = mkseries(T2, 2000, 0, 0.1);
endif 

## trying to do county population
if 0
    lb_spec_p = zeros(36);
    ub_spec_p = ones(36)+2;
    eq_spec_p = blkdiag(diag(nan(17,1), -1), diag(nan(17,1), -1));
    eq_spec_p(1, 20:26) = nan;
    eq_spec_p(19, 20:26) = nan;
    P=Pop();
    p1234 = P.pop{1234};
endif


if 0

    #########################
    ## caswell code -- works just like the book!
    z = Pc(:,2:end)(:);
    nonzero = [1 2 5 6 7 9];
    M = [];
    for ii=1:5
        N = kron(Pc(:,ii)', eye(3));
        m = N(:, nonzero);
        M = [M; m];
    endfor
    G = M'*M; 
    f = -M'*z;
    C = diag(-ones(6,1));
    C(end+1,:) = [1 1 0 0 0 0];
    C(end+1,:) = [0 0 1 1 0 0];
    C(end+1,:) = [0 0 0 0 0 1];
    b = [0 0 0 0 0 0 1 1 1]';
    [phat cas.obj cas.info cas.lambdas] = qp([], G, f, [], [], [], [], -inf(9,1), C, b);
    a = zeros(9,1);             # empty matrix to be filled
    a(nonzero) = phat;
    a = reshape(a, 3, []); 
    cas._ssq       = sumsq(Pc(:,2:end)(:));
    cas._term1     = f'*phat;
    cas._term2     = 0.5*phat'*G*phat;
    cas.ssqmod     = sumsq((a*Pc(:,1:end-1))(:));

    ########################
    ## sqp, full matrix, and regularized...
    kI = eye(rows(Pc));
    M  = kron(Pc(:,1:end-1)', kI);        
    z  = Pc(:,2:end)(:);
    a_stuff = struct();
    a_stuff.H = (M'*M);                 # note: our function multiplies H by 0.5 so we don't here
    a_stuff.q = -(M'*z);

    _mysqp = struct();

    ## equality -- sparsity in this case
    x = mkconeq(eq_spec);
    a_stuff.A_eqm = x.m;
    a_stuff.b_eqv = x.v;

    ## inequality, incomplete currently XXX (no above positive by element)
    x = mkconineq(ineq_spec, 1); 
    a_stuff.A_in = x.m;
    a_stuff.A_lb = x.v;

    a_stuff.lbvec = zeros(9,1);

    a_stuff.x0   = rand(9,1) - 0.5; 

    [_mysqp.m _mysqp.D] = wmqp('sqp1', a_stuff);
    _mysqp.m = reshape(_mysqp.m, 3, 3);

endif


#{
[m D] = fitm(p(:,1:end-1), p(:,2:end), [], [], [], [], [], 'treg', 0.0, 'median', 0);
merr(end+1) = sumsq(p(:,end) - m*p(:, end-1));
mobj(end+1) = D.obj;

[m D] = fitm(p(:,1:end-1), p(:,2:end), [], [], [], [], [], 'treg', 0.0, 'median', 1);
merr(end+1) = sumsq(p(:,end) - m*p(:, end-1));
mobj(end+1) = D.obj;

[m D] = fitm(p(:,1:end-1), p(:,2:end), eq_spec1, [], [], [], [], 'treg', 0.0, 'median', 1);
merr(end+1) = sumsq(p(:,end) - m*p(:, end-1));
mobj(end+1) = D.obj;

[m D] = fitm(p(:,1:end-1), p(:,2:end), [], {lc_spec1}, [], lb_spec1, ub_spec1, 'median', 1);
merr(end+1) = sumsq(p(:,end) - m*p(:, end-1));
mobj(end+1) = D.obj;

[m D] = fitm(p(:,1:end-1), p(:,2:end), [], [], [], lb_spec1, ub_spec1, 'median', 1);
merr(end+1) = sumsq(p(:,end) - m*p(:, end-1));
mobj(end+1) = D.obj;
  
[m D] = fitm(p(:,1:end-1), p(:,2:end), [], {lc_spec1}, [], lb_spec1, ub_spec1, 'median', 1);
merr(end+1) = sumsq(p(:,end) - m*p(:, end-1));
mobj(end+1) = D.obj;

[m D] = fitm(p(:,1:end-1), p(:,2:end), eq_spec1, {lc_spec1}, [], lb_spec1, ub_spec1, 'median', 1);
merr(end+1) = sumsq(p(:,end) - m*p(:, end-1));
mobj(end+1) = D.obj;

[m D] = fitm(p(:,1:end-1), p(:,2:end), [], {lc_spec1}, [], lb_spec1, ub_spec1, 'treg', 100, 'median', 1);
merr(end+1) = sumsq(p(:,end) - m*p(:, end-1));
mobj(end+1) = D.obj;

T = [0 0 0 1; 1 0 0 0; 0 1 0 0; 0 0 1 0];
T = [0 0.7 0.7 0; .95 0 0 0; 0 .95 0 0; 0 0 .5 0];


#}
