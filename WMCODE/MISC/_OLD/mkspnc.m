function out = mkspnc(indata, kwig=0.0, neigs=[])
    ## Creates constraints from input data. "A" and "b" enforce a
    ## linear combination of principal components derived from the
    ## SVD. "A_in", "a_lb", "a_ub" enforce inequality bounds on the
    ## component multipliers.
    ## 
    ## kwig defines multiplier on upper and lower bounds
    ## 
    ## neigs defines the number of components returned

    ## prep data
    _c = columns(indata);
    _r = rows(indata);
    if isempty(neigs) 
        neigs = min(_r, _c);
    endif
        
    ## Eigen wiggle clean and warn
    if kwig ~= 0.0
        kwig = abs(kwig);
        warning("Using non zero k-wiggle: %f", kwig);
    endif
    
    ## SVD magic with centering -- 
    davg = mean(indata);
    [u s v] = svd(indata - repmat(davg, rows(indata), 1));    
    
    ## fix negatives
    if all(u<=0)(1)
        u = -1*u; 
    endif
    
    ## keep only the good ones
    u_ = u(:,1:neigs);
    s_ = s(1:neigs,1:neigs);
    v_ = v'(1:neigs,:);

    ## scale S if using neigs
    if neigs < columns(indata)
        out.scs = cumsum(diag(s)/sum(diag(s)));
        sf =  1/(out.scs(neigs));
        s_ = s_ * sf;
    endif
    
    ## ... build eq coeff matrices
    coeffs = (s_ * v_)';
    out.A = [eye(_c), coeffs];
    out.b = davg';

    ## ... ineq matrices, from u using k-wiggle
    kmin = min(u_)';
    kmax = max(u_)';
    kspr = abs(kmin - kmax) * kwig;
    out.a_lb = kmin - kspr;
    out.a_ub = kmax + kspr;
    out.A_in = [zeros(neigs, _c), eye(neigs)];

    ## H for qp() later
    out.H = eye(columns(out.A));
    out.H((_c+1):end, (_c+1):end) = 0;

    ## more info is better...
    out.u = u; out.s = s; out.v = v; out._c = _c; out._r = _r; out.u_ = u_;

endfunction    


#[x, obj, info, lambda] = qp (x0, h, q, a, b, lb, ub, a_lb, a_in, a_ub)

%!demo
%!   "Fertility"
%!   ages1 = 15:44;
%!   ages5 = 15:5:40;
%!   F = Fert;
%!   neigs = 4;
%!   L_i = neigs*2 - 1;
%!   Cs = mkspnc(F.Fx5(find(F.r_year >= 1950),:), 0, neigs);
%!   disp(Cs.scs(1:neigs)');

%!   for fi = 1:rows(F.Fx5)
%!       target = [F.Fx5(fi,:)'; zeros(neigs,1)];
%!       [X, OBJ, INFO, LAMBDA] = qp([], Cs.H, -1*target, Cs.A, Cs.b, [], [], Cs.a_lb, Cs.A_in, Cs.a_ub);
%!       disp(INFO);
%!       plot(ages5, F.Fx5(fi,:)', 'b');
%!       ylim([0 1.2]);
%!       hold on;
%!       plot(ages5, X(1:end-neigs), 'r');
%!       hold off;
%!       disp(LAMBDA(end-L_i:end)');
%!       printf ("i: %i, tfr: %f\n", fi, sum(F.Fx5(fi,:)));
%!       input ("press ","s");
%!   endfor


%!demo
%!   "Survival"
%!   ages5 = 5:5:85;
%!   S = Surv();
%!   neigs = 4;
%!   L_i = neigs*2 - 1;
%!   Cs = mkspnc(S.Surv(:,:), 0, neigs);
%!   disp(Cs.scs(1:neigs)');

%!   for fi = 1:rows(S)
%!       target = [S.Surv(fi,:)'; zeros(neigs,1)];
%!       [X, OBJ, INFO, LAMBDA] = qp([], Cs.H, -1*target, Cs.A, Cs.b, [], [], Cs.a_lb, Cs.A_in, Cs.a_ub);
%!       disp(INFO);
%!       subplot(2,1,1);
%!       plot(ages5, S.Surv(fi,1:17)', 'b');
%!       ylim([0.3 1.10]);
%!       xlim([0 90]);
%!       line ([0 85], [0 0 ]);
%!       hold on;
%!       plot(ages5, X(1:17), 'r');
%!       hold off;
%!       subplot(2,1,2);
%!       plot(ages5, S.Surv(fi,18:34)', 'b');
%!       ylim([0.3 1.10]);
%!       xlim([0 90]);
%!       line ([0 85], [0 0 ]);
%!       hold on;
%!       plot(ages5, X(18:34), 'r');
%!       hold off;
%!       disp(LAMBDA(35:end)');
%!       printf ("i: %i, yr: %i, geo: %s\n", fi, S.r_year(fi), S.r_geo{fi});
%!       input ("press ","s");
%!   endfor



%!demo
%!   "WA only migration, using only model counties"
%!   ages5 = 0:5:85;
%!   MM = Mig();
%!   neigs = rows(MM.WA90_model);
%!   L_i = neigs*2 - 1;
%!   Cs = mkspnc(MM.WA90_model(:,:), 100);

%!   for fi = 1:rows(MM.WA90)
%!       target = [MM.WA90(fi,:)'; zeros(neigs,1)];
%!       [X, OBJ, INFO, LAMBDA] = qp([], Cs.H, -1*target, Cs.A, Cs.b, [], [], Cs.a_lb, Cs.A_in, Cs.a_ub);
%!       disp(INFO);
%!       plot(ages5, MM.WA90(fi,:)', 'b');
%!       grid on;
%!       ylim([-0.6 1.0]);
%!       line ([0 90], [0 0 ]);
%!       hold on;
%!       plot(ages5, X(1:end-neigs), 'r');
%!       hold off;
%!       disp(LAMBDA(end-L_i:end)');
%!       printf ("i: %i\n", fi);
%!       input ("press ","s");
%!   endfor



%!demo
%!   "National migration"
%!   ages5 = 0:5:85;
%!   popmin = 2500;
%!   MM = Mig();
%!   popmin_i = find(sum(MM.mfpop,2) >= popmin);
%!   neigs = 12;
%!   Cs = mkspnc(MM.NetMig(popmin_i,:), 0, neigs);
%!   ORi = strmatch('41001', MM.r_geo):strmatch('41071', MM.r_geo);
%!   WAi = strmatch('53001', MM.r_geo):strmatch('53077', MM.r_geo);
%!   CAi = strmatch('06001', MM.r_geo):strmatch('06115', MM.r_geo);


%!   for fi = WAi
%!       target = [MM.NetMig(fi,:)'; zeros(neigs,1)];
%!       [X, OBJ, INFO, LAMBDA] = qp([], Cs.H, -1*target, Cs.A, Cs.b, [], [], Cs.a_lb, Cs.A_in, Cs.a_ub);
%!       disp(INFO);

%!       subplot(2,1,1);
%!       plot(ages5, MM.NetMig(fi,1:18)', 'b');
%!       grid on;
%!       ylim([-0.6 1.0]);
%!       line ([0 90], [0 0]);
%!       hold on; plot(ages5, X(1:18), 'r', 'linewidth', 1.75); hold off;

%!       subplot(2,1,2);
%!       plot(ages5, MM.NetMig(fi,19:36)', 'b');
%!       grid on;
%!       ylim([-0.6 1.0]);
%!       line ([0 90], [0 0]);
%!       hold on; plot(ages5, X(19:36), 'r', 'linewidth', 1.75); hold off;

%!       disp(LAMBDA(37:end)');
%!       printf ("i: %i, pop: %i, fips: %5.5s\n", fi, sum(MM.mfpop(fi,:)), MM.r_geo{fi,1});
%!       input ("press ","s");
%!   endfor



