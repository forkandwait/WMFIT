function [c stuff] =  mkconeq(ccell, conl=[])
    ## make per-cell constraints for qp.
    ##
    ## ccell creates simple per-cell equality constraints
    ##
    ## conl creates multi-cell constraints by specifying a matrix of
    ## coefficients and a scalar to sum 
    ## 
    ## {{coeff1, sum1}, {coeff2, sum2}} 



    ## For between constrains, mkconbw 
      
    ccell = ifelse(isnumeric(ccell), num2cell(ccell), ccell);

    ## set up the output matrices
    matl = numel(ccell);
    m = nan(0, matl);
    v = nan(0, 1);

    ## create the per-cell equalities
    ii = 1;
    for _cc = ccell(:)'
        cc = _cc{1};
        if isempty(cc) || isnan(cc)
           ;
        else
            cc = _cc{1};
            crow = zeros(1, matl);
            crow(ii) = 1;
            m(end+1, :) = crow;
            v(end+1, 1) = cc;
        endif
        ii++;    
    endfor

    ## create the multi-cell equalities
    conl = conl(:);
    for ii = 1:length(conl)
        _cc = conl{ii};
        _m = _cc{1}(:);
        if numel(_m) ~= matl
            error("mkcon: weird size for row sum vector: should have numel == numel of transition matrix. %i ~= %i.", 
                  numel(_m), numel(ccell));
        endif
        m(end+1, :) = _m;
        v(end+1, 1) = _cc{2};
    endfor

    ## package up and return
    c = struct();
    c.m = m;
    c.v = v;
endfunction

%!test
%!  c = mkcon ({ 1, {}, 5, {[1 1; 1 0], 3} });
%!  assert(c.m, [1 0 0 0; 0 0 1 0; 1 1 1 1]);
%!  assert(c.v, [1; 5; 3]);

%!test
%!  ceq = mkcon ({0 {} .5 0; {} 0 0 0; 0 {} 0 0; 0 0 {} 0});  

%!test
%!  ceq = mkcon ({0 {[0 0 1 0; 0 0 0 0; 0 0 0 0; 0 0 0 0], 0} {}  0; {} 0 0 0; 0 {} 0 0; 0 0 {} 0});  


%!#test -- this doesn't test for fail like it should
%!  fail(mkcon ({ 1, {}, 5, {[1 2; 3 0], [6 3]} }), "mkcon");
