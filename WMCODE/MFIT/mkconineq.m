function [c stuff] =  mkconineq(ccell, issingle=0)
    ## constrain lincombos to be inequality between numbers.  Spec
    ##     format: { {matrix, lb, ub}, {...}, }
    ## 
    ## Redundant with lb and ub constraints, but still useful, I think
         
    mf = mfilename;
    ccell = ccell(:);

    matl = numel(ccell{1}{1});
    m    = nan(0, matl);
    vl   = nan(0, 1);
    v2   = nan(0, 1);

    ## do the matrix building 
    for ii = 1:length(ccell);
        cc = ccell{ii};
        if isempty(cc)
            ;
        elseif (iscell(cc) && all(size(cc) == [1 3]) && numel(cc{2}) == 1 && numel(cc{3}) == 1)
            _m = cc{1}(:);
            if numel(_m) ~= matl
               error("%s: weird size for row sum vector: should have numel == numel of transition matrix. %i ~= %i.", 
                     mf, numel(_m), numel(ccell));
            endif
            m(end+1, :)  = _m;
            v1(end+1, 1) = cc{2};
            v2(end+1, 1) = cc{3};
            if cc{2} >= cc{3} 
               error("%s: lower bound >= upper bound %f %f", mf, cc{2}, cc{3});
            endif
        else
            error("%s: cell must be special 3 element cell array", mf);
        endif
    endfor 

    ## assemble the output struct
    c = struct();
    if issingle                 # combined lb and ub into lb, using negatives
        c.m = vertcat(m,  -1*m);
        c.v = vertcat(v1, -1*v2);
    else 
        c.m = m;
        c.v1 = v1;
        c.v2 = v2;
    endif

endfunction

%!test
%!  ineqspec = {
%!        {[1 1; 0 0],    0, 1},
%!        {[1 1; 1 0], -Inf, 3}
%!             };
%!  c = mkconineq (ineqspec);

%!test
%!  ineqspec = { {[1 1 0 0], 0, 1} };
%!  c = mkconineq (ineqspec);


%!#test -- this doesn't test for fail like it should
%!  fail(mkconineq ({ 1, {}, 5, {[1 2; 3 0], [6 3]} }), "mkcon");
