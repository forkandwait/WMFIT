## experiment with sqp regularized stuff

function out = show4(m, d=4)
    out = floor(m*10^4)*10^-4;
endfunction

## simple, based on caswell
if 1

    ## basic data
    T  = [0.0000 0.0000 0.0143; 
          1.0000 0.4611 0.0000; 
          0.0000 0.5389 0.9664];
    rand("seed", 1);
    p  = mkseries(T, 5, [], [], [5.32 24.84 115.5]');
    rand("seed", 1);
    p2 = mkseries(T, 50, [], [], [5.32 24.84 115.5]');
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
    [T0_fit V0 D0] = fitm(p2(:,1:end-1), p2(:,2:end), eq_spec, [], lb_spec, [], ineq_spec, 'algo', 'qp');
    [T1_fit V1 D1] = fitm(p2(:,1:end-1), p2(:,2:end), eq_spec, [], lb_spec, [], ineq_spec, 'algo', 'sqp1');
    [T2_fit V2 D2] = fitm(p2(:,1:end-1), p2(:,2:end), [], [], [], [], [], 'algo', 'sqp1');
    [T3_fit V3 D3] = fitm(p2(:,1:end-1), p2(:,2:end), [], [], [], [], [], 'algo', 'sqp2', 'alpha', 1);
    [T4_fit V4 D4] = fitm(p2(:,1:end-1), p2(:,2:end), [], [], [], [], [], 'algo', 'sqp2', 'alpha', 0.25);
endif

## with linear combination
if 0
    T = [ 
          0.00000   0.63450   0.63450   0.00000
          0.90000   0.00000   0.00000   0.00000
          0.00000   0.75000   0.00000   0.00000
          0.00000   0.00000   0.50000   0.00000
    ];
    
    eq_spec1 = zeros(4);
    eq_spec1(1, 2:3) = nan;
    eq_spec1(2,1) = eq_spec1(3,2) = eq_spec1(4,3) = nan;
    
    lc_spec{1}      = struct();
    lc_spec{1}.name = 'fert';
    lc_spec{1}.idx  = [1 2; 1 3]; 
    lc_spec{1}.vset = { struct('lc', [.5 .5], 'vlb', .5, 'vub', 4) };
    lc_spec{1}.b    = [0 0];
    
    lb_spec1 = zeros(4);
    ub_spec1 = ones(4);
    
    lb_spec2 = repmat(-1, 4, 4);
    ub_spec2 = repmat(1, 4, 4);

    p = mkseries(T, 2000, 0, 0.1);

    [T1_fit D1] = fitm(p(:,1:end-1), p(:,2:end), eq_spec1, [], lb_spec1, [], [], 'algo', 'sqp2', 'alpha', 1.0);
    [T4_fit D4] = fitm(p(:,1:end-1), p(:,2:end), eq_spec1, [], lb_spec1, [], [], 'algo', 'sqp2', 'alpha', 0.5);


endif 


## with human population 
if 0
   if !exist('P')
       P = Pop();
       spoki = strmatch('53063', P.r_geo);
       p = P.pop{spoki};
   endif

   ## define the sparsity of human populations:
   _foo = diag(ones(1,17), -1) + eye(18);
   eq_spec = blkdiag(_foo, _foo);
   eq_spec(1, [19 21:28]) = 1;
   eq_spec(1, [19 21:28]) = 1;
   eq_spec(eq_spec==1)    = nan;

   disp("    T3");
   [T3_fit D3] = fitm(p(:,1:end-1), p(:,2:end), eq_spec, [], [], [], [], 'algo', 'sqp2', 'alpha', 1.00);

   #{

   disp("    T4");
   [T4_fit D4] = fitm(p(:,1:end-1), p(:,2:end), [], [], [], [], [], 'algo', 'sqp2', 'alpha', 0.99);

   disp("    T1");
   [T1_fit D1] = fitm(p(:,1:end-1), p(:,2:end), [], [], [], [], [], 'algo', 'sqp2', 'alpha', 1.0);

   #}
   
endif
