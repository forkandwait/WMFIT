function constr = mkconsvd (sdata, fdata, mdata, 
                            sneigs=4, fneigs=4, mneigs=12, 
                            swig=0, fwig=100, mwig=100)
    ## make constraints from data and svd.  Modify G and f input for
    ## qp() to deal with extra floating variables.

    constr = struct();

    ## get svd constraints
    constr._scon = mkspnc(sdata, swig, sneigs);    
    constr._fcon = mkspnc(fdata, fwig, fneigs);    
    constr._mcon = mkspnc(mdata, mwig, mneigs);

	## get matrix infrastructure
	mstuff = mkmat("twosex");

    ## XXXX  Won't be using pure bounds -- do everything in bounds matrices
    constr.lb = [];
    constr.ub = [];

    #### Matrix infrastructure ################
    ## set up some indices for the matrix craziness
    _fvarcnt = sneigs + fneigs + mneigs; # number of "free" variables to control demography
    _mvarcnt = length(mstuff.nz); # number of "matrix" variables in output leslie matrix

    ## ... finding components in the constraints matrices -- finish
    ## with three vectors of numbers from find()
    _sidx = zeros(1, _mvarcnt + _fvarcnt);
    _fidx = zeros(1, _mvarcnt + _fvarcnt);
    _midx = zeros(1, _mvarcnt + _fvarcnt);
    _sidx(1, (_mvarcnt+1):(_mvarcnt+sneigs)) = 1;
    _fidx(1, (_mvarcnt+sneigs+1):(_mvarcnt+sneigs+fneigs)) = 1;
    _midx(1, (_mvarcnt+sneigs+fneigs+1):(_mvarcnt+_fvarcnt)) = 1;
    _sidx = find(_sidx);
    _fidx = find(_fidx);
    _midx = find(_midx);

    ## ... finding components in the leslie matrix
    _msidx = [mstuff._sim, mstuff._sif]; # survival 
    _mfidx = [mstuff._fim, mstuff._fif]; # fertility
    _mmidx = [mstuff._dim, mstuff._dif]; # migration

	#### Equality constraints. ################

    ## .. empty setup to be added to ..
	constr.Aeq = zeros(0, _mvarcnt + _fvarcnt);
	constr.beq = zeros(0, 1);

    ##     .. survival ..
	_Aeq = zeros(0, _mvarcnt + _fvarcnt);
	_beq = zeros(0, 1);
    for _ri = 1:rows(constr._scon.A)
        _Aeq(_ri, _msidx(_ri)) = 1;
        _Aeq(_ri, _sidx)       = constr._scon.A(_ri, end-sneigs+1:end);
        _beq(_ri,1)            = constr._scon.b(_ri);
    endfor
    constr.Aeq = vertcat(constr.Aeq, _Aeq);
    constr.beq = vertcat(constr.beq, _beq);

    ##     .. fertility ..
    warning("split fertility .4886!!!");
	_Aeq = zeros(0, _mvarcnt + _fvarcnt);
	_beq = zeros(0, 1);
    for _ri = 1:rows(constr._fcon.A)
        _Aeq(_ri, _mfidx(_ri)) = 1;
        _Aeq(_ri, _fidx)       = constr._fcon.A(_ri, end-fneigs+1:end);
        _beq(_ri,1)            = constr._fcon.b(_ri);
    endfor
    constr.Aeq = vertcat(constr.Aeq, _Aeq);
    constr.beq = vertcat(constr.beq, _beq);

    ##     .. migration ..
    if 1
	    _Aeq = zeros(0, _mvarcnt + _fvarcnt);
	    _beq = zeros(0, 1);
        for _ri = 1:rows(constr._mcon.A)
            _Aeq(_ri, _mmidx(_ri)) = 1;
            _Aeq(_ri, _midx)       = constr._mcon.A(_ri, end-mneigs+1:end);
            _beq(_ri,1)            = constr._mcon.b(_ri);
        endfor
        constr.Aeq = vertcat(constr.Aeq, _Aeq);
        constr.beq = vertcat(constr.beq, _beq);
    endif

    ##### Inequality constraints ..  ################
    constr.A_in = zeros(0, _mvarcnt + _fvarcnt);
    constr.a_ub = zeros(0,1);
    constr.a_lb = zeros(0,1);

    ## ... survival ..
    _A_in = zeros(0, _mvarcnt + _fvarcnt);
    _a_ub = zeros(0,1);
    _a_lb = zeros(0,1);
    for _ri = 1:rows(constr._scon.A_in)
        _A_in(_ri, _msidx(_ri)) = 1;
        _A_in(_ri, _sidx)       = constr._scon.A_in(_ri, end-sneigs+1:end);
        _a_lb(_ri,1)            = constr._scon.a_lb(_ri);
        _a_ub(_ri,1)            = constr._scon.a_ub(_ri);
    endfor
    constr.A_in = vertcat(constr.A_in, _A_in);
    constr.a_lb = vertcat(constr.a_lb, _a_lb);
    constr.a_ub = vertcat(constr.a_ub, _a_ub);

    ## ... fertility ..
    _A_in = zeros(0, _mvarcnt + _fvarcnt);
    _a_ub = zeros(0,1);
    _a_lb = zeros(0,1);
    for _ri = 1:(rows(constr._fcon.A_in))
        _A_in(_ri, _mfidx(_ri)) = 1;
        _A_in(_ri, _fidx)       = constr._fcon.A_in(_ri, end-fneigs+1:end);
        _a_lb(_ri,1)            = constr._fcon.a_lb(_ri);
        _a_ub(_ri,1)            = constr._fcon.a_ub(_ri);
    endfor
    constr.A_in = vertcat(constr.A_in, _A_in);
    constr.a_lb = vertcat(constr.a_lb, _a_lb);
    constr.a_ub = vertcat(constr.a_ub, _a_ub);

    ## ... migration ..
    if 1
        _A_in = zeros(0, _mvarcnt + _fvarcnt);
        _a_ub = zeros(0,1);
        _a_lb = zeros(0,1);
        for _ri = 1:rows(constr._mcon.A_in)
            _A_in(_ri, _mmidx(_ri)) = 1;
            _A_in(_ri, _midx)       = constr._mcon.A_in(_ri, end-mneigs+1:end);
            _a_lb(_ri,1)            = constr._mcon.a_lb(_ri);
            _a_ub(_ri,1)            = constr._mcon.a_ub(_ri);
        endfor
        constr.A_in = vertcat(constr.A_in, _A_in);
        constr.a_lb = vertcat(constr.a_lb, _a_lb);
        constr.a_ub = vertcat(constr.a_ub, _a_ub);
    endif

endfunction

%!demo
%!    fipsi = randi([1 3120])
%!    S = Surv();
%!    F = Fert();
%!    M = Mig();
%!    P = Pop();
%!    constr = mkconsvd(S.Surv(250:256,:), F.Fx5_2(248:255,:), M.NetMig); 
%!    [A, obj, info, lambda] = wmh(P.pop{fipsi}, constr);
%!    disp([round(A*P.pop{fipsi}(:,end-1)),  P.pop{fipsi}(:,end)]);
%!    plot(round(A*P.pop{fipsi}(:,end-1))); hold on; plot(P.pop{fipsi}(:,end), 'r'); hold off;


%!demo
%!    S = Surv();
%!    F = Fert();
%!    M = Mig();
%!    P = Pop();
%!    constr = mkconsvd(S.Surv(250:256,:), F.Fx5_2(248:255,:), M.NetMig); 
%!    for fi = 2932:2970
%!        [A, obj, info, lambda] = wmh(P.pop{fi}(:,1:(end-2)), constr);
%!        disp([round(A*P.pop{fi}(:,end-1)),  P.pop{fi}(:,end)]);
%!        figure(1);
%!        plot(round(A*P.pop{fi}(:,end-1)), 'r'); hold on; plot(P.pop{fi}(:,end-1), 'g'); plot(P.pop{fi}(:,end), 'b'); hold off; grid on;
%!        figure(2);
%!        plot(round(A^2*P.pop{fi}(:,end-2)), 'r'); hold on; plot(P.pop{fi}(:,end-2), 'g'); plot(P.pop{fi}(:,end), 'b'); hold off; grid on;
%!        input(sprintf("%s, %i\n", P.fips{fi}, info.solveiter));
%!    endfor


%!demo
%!    fipsi = 2969;
%!    cfipsi = [2960 2979];
%!    S = Surv();
%!    F = Fert();
%!    M = Mig();
%!    cmig =  M.NetMig(cfipsi,:); cmig(:,1) = 0;
%!    P = Pop();
%!    constr = mkconsvd(S.Surv(250:256,:), F.Fx5_2(248:255,:), cmig, 4, 4, 2, 1, 1, 1); 
%!    [A, obj, info, lambda] = wmh(P.pop{fipsi}(:,1:(end-2)), constr);
%!    disp([round(A*P.pop{fipsi}(:,end-1)),  P.pop{fipsi}(:,end)]);
%!    plot(round(A*P.pop{fipsi}(:,end-1)), 'r'); hold on; plot(P.pop{fipsi}(:,end), 'b'); hold off;
