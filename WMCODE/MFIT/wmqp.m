function [x D] = wmqp(spec='qp', argstr)
    ## execute quadratic programming for woods method, different types
    ##     for different id's.
    ##
    ## argstr is a struct to hold all the various components:
    ##
    ##    x0, H, q, A_eqm, b_eqv, lbvec, ubvec, A_lb, A_in, A_ub
    ## 
    ## x    -- the optimized values
    ##
    ## D    -- debug dump of all variables, including all the
    ##             optimization information
         
    argstr = mkargstr(argstr);

    switch spec
        case 'qp'
            [x obj info lambda] = qp(argstr.x0, argstr.H, argstr.q, 
                                     argstr.A_eqm, argstr.b_eqv, 
                                     argstr.lbvec, argstr.ubvec, 
                                     argstr.A_lb, argstr.A_in, argstr.A_ub);
            if info.info ~= 0
                warning("wmqp: bad info for qp: %i", info.info);
            endif 
        case 'sqp1'
             phi = mkphi1(argstr.H, argstr.q); 
             g   = ifelse(isempty(argstr.A_eqm), [], mkg1(argstr.A_eqm, argstr.b_eqv));
             h   = ifelse(isempty(argstr.A_in),  [], mkh1(argstr.A_in, argstr.A_lb));
             [x, obj, info, iter, nf, lambda] = sqp(argstr.x0, phi, g,   h, argstr.lbvec, argstr.ubvec); # real thing
             if info ~= 101
                warning("wmqp: bad info for sqp: %i. 102=\"BGFS update failed\", 103=\"Max iterations\"", info);
             endif
        case 'sqp2'
             phi = mkphi2(argstr.H, argstr.q, argstr.alpha); 
             g   = ifelse(isempty(argstr.A_eqm), [], mkg1(argstr.A_eqm, argstr.b_eqv));
             h   = ifelse(isempty(argstr.A_in),  [], mkh1(argstr.A_in, argstr.A_lb));
             [x, obj, info, iter, nf, lambda] = sqp(argstr.x0, phi, g, h, argstr.lbvec, argstr.ubvec);
             if info ~= 101
                warning("wmqp: bad info for sqp: %i. 102=\"BGFS update failed\", 103=\"Max iterations\"", info);
             endif
        otherwise
            error("wmqp: unknown spec: %s", spec);
    endswitch

    ## calculate norm from x0

    ## return all the intermediate stuff if 2 output args
    if nargout == 2
       D = allv();
    endif

endfunction


## objective
function phi = mkphi1(H, q)
    phi{1} = @(x)(0.5*x(:)'*H*x(:) + q(:)'*x(:));
    phi{2} = @(x)(H*x(:) + q(:));
    phi{3} = @(x)(H);
endfunction

## objective, try to l-2 regularize...
function phi = mkphi2(H, q, alpha)
    if alpha < 0 || alpha > 1
        error("wmqp: alpha must be between 0 and 1. Is: %f.", alpha);
    endif
    phi{1} = @(x)(alpha*(0.5*x(:)'*H*x(:) + q(:)'*x(:)) + (1-alpha)*sum(abs(x(:))));
    phi{2} = @(x)(alpha*(H*x(:) + q(:)) + (1-alpha)*sign(x(:)));
    phi{3} = @(x)(H);
endfunction

## g(x) = 0
function g   = mkg1(A_eqm, b_eqv)
    g{1}   = @(x)(A_eqm * x - b_eqv);
    g{2}   = @(x)(A_eqm); 
endfunction

## h(x) >= 0
function h   = mkh1(A_in, A_lb)
    h{1}   = @(x)(A_in * x - A_lb);
    h{2}   = @(x)(A_in); 
endfunction

function out = mkargstr (argstr)
    out = argstr;

    ## populate argument structure with default []'s
    for _s = {'x0', 'H', 'q', 'A_eqm', 'b_eqv', 'lbvec', 'ubvec', 'A_lb', 'A_in', 'A_ub', 'alpha', 'opstr'}
        if !isfield(out, _s{1})
            out = setfield(out, _s{1}, []);
        endif
    endfor
    
    ## deal with empty lbvec/ ubvec
    if isempty(out.lbvec) && !isempty(out.ubvec) 
       warning ("wmqp: compensating for empty ubvec");
       out.lbvec = inf(size(argstr.ubvec));
    elseif !isempty(out.lbvec) && isempty(out.ubvec) 
       warning ("wmqp: compensating for empty lbvec");
       out.ubvec = inf(size(argstr.lbvec));
    endif

    ## deal with empty alpha -- set to 1,  for all fit
    if isempty(out.alpha) 
       out.alpha = 1.0;
       warning("wmqp: no alpha specified, setting to zero");
    endif

endfunction


%!test 
%!    a       = struct();
%!    a.x0    = rand(3,1);
%!    a.H     = eye(3);
%!    a.q     = [0 0 0]';
%!    a.A_eqm = [1 0 0; 1 1 0];
%!    a.b_eqv = [10, 1000]';
%!    a.A_in  = [1 1 1];
%!    a.A_lb  = 4000;
%!    [x D] = wmqp(spec='sqp1', a);
%!    disp(x);

%!test 
%!    a       = struct();
%!    a.x0    = rand(3,1);
%!    a.H     = eye(3);
%!    a.q     = [0 0 0]';
%!    a.A_eqm = [1 0 0; 1 1 0];
%!    a.b_eqv = [10, 1000]';
%!    a.A_in  = [];
%!    a.A_lb  = [];
%!    a.alpha = .5;
%!    [x D] = wmqp(spec='sqp2', a);
%!    disp(x);

