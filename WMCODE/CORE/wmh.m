function [A, stuff] = wmh(inpop, constr=[], varargin)
    ## [A, stuff] = wmh(inpop, constr=[], named_params)

    ## Double check file, as this function is a  work in progress and
    ## parameters change regularly

    ## named options
    [reg prop] = parseparams(varargin);
    Opts = struct(prop{:});

    ## struct to hold random useful things
    debug = struct();

    ## get matrix infrastructure
    mat = mkmat("twosex");

    ## get constraints		
    if isempty(constr) 
       cratio = mean(inpop(5,:) ./ inpop(6,:));
       if cratio >= 1.75
           ## warning("%s: using college constraints. cratio = %f.", mfilename(), cratio);
           constr = mkconwm("college");
       else
           constr = mkconwm();
       endif
    endif

    ## vec() applied to noutx
    popv = (inpop(:, 2:end))(:);

    ## generate matrix M
    M = [];
    for i = 1:(columns(inpop)-1)
	N = kron(inpop(:,i)', eye(rows(inpop)));
	m = N(:, mat.nz);
	M = [M; m];
    endfor

    ## matrix G and f
    G = (M' * M);
    f = -M' * popv;

    ## fix size of G and f if constraints.
    _fvarcnt = columns(constr.Aeq) - columns(G);
    G = blkdiag(G, zeros(_fvarcnt));
    f = vertcat(f, zeros(_fvarcnt, 1));

    ## Do quadratic optimization -- three times with sleep ;)
    [phat, stuff.obj, stuff.info, stuff.lambda] = qp([], G, f, constr.Aeq, constr.beq,
                                                     constr.lb, constr.ub, 
                                                     constr.a_lb, constr.A_in, constr.a_ub);
    ii = 1;
    while (any(isnan(phat)) || (stuff.info.info == 3)) 
	    if ii >= 4
	        warning("%s: not trying anymore", mfilename);
	        break;
	    end
	    warning("%s: try again -- time %i, sleep: %i.", mfilename, ii, 2^ii);
	    sleep(2^ii);
	    [phat, stuff.obj, stuff.info, stuff.lambda] = qp([], G, f, constr.Aeq, constr.beq,
                                                         constr.lb, constr.ub, 
                                                         constr.a_lb, constr.A_in, constr.a_ub);
	    ii += 1;
    end

    ## warn on non-convergent
    if stuff.info.info ~= 0
        warning("wmh: Problem with qp. info: %i, solveiter: %i", stuff.info.info, stuff.info.solveiter);
    end

    ## new projection matrix.
    A = zeros(rows(inpop));
    A(mat.nz) = phat(1:length(mat.nz)); # assign only cell variables, leave free coordination variables out

    ## post-fill matrix with stuff from varargin
    if isfield(Opts, "Fert")
       _fert = getfield(Opts, "Fert");
       if ! all(size(_fert) == [9 1])
          error("wmh: Fert should be 9 x 1, is %i x %i.", size(_fert));
       endif
       _fert = [_fert * (1-0.4886), _fert * 0.4886];
       A(1,1) = _fert(1,1);
       A(19,19) = _fert(1,2);
       A(1, 21:28) = _fert(2:9,1);
       A(19,21:28) = _fert(2:9,2);
    endif
       

endfunction
