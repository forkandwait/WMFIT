function [T DEBUG] = mksmth(smth_spec)
    ## mksmth returns a 2nd derivative regularization matrix (Boyd)
    ## 
    ## smth_spec is a matrix with unique non-zero integers defining
    ## sections of a matrix to try to smooth together
    ## 
    ## [0 1 1 0;
    ##  2 0 0 0;
    ##  0 2 0 0;
    ##  0 0 2 0];
    ## 
    ## This matrix separately smooths the top row and the
    ## subdiagonal.
    ## 
    ## The implementation is provisional, at best.  The most obvious
    ## problem is with the order of the sub-spec and how it defines
    ## the boundaries. It is arbitrarily in the same order as the
    ## vectorized matrix (row major), with the first and last of the
    ## given spec number defining the boundary of the smoothing range;
    ## each middle entry is [-1 2 -1], but the the boundaries are [0 2
    ## -1]. XXX verify, rework, etc.
         
    error("mksmth: not ready for prime time -- need sqp...");
         
    ## make and fill output matrix get and treat the various unique indices
    tsz = size(smth_spec);
    T = zeros(prod(tsz));
    items = unique(smth_spec(:));

    for ii = items(2:end)' 
        celli = find(smth_spec(:)==ii);
        stp = length(celli);
        for jj = 1:stp
            T(celli(jj), celli(jj))   = -2;
            switch jj 
                case 1 
                    T(celli(jj), celli(jj)+1) = -1;
                case stp
                    T(celli(jj), celli(jj)-1) = -1;
                otherwise
                    T(celli(jj), celli(jj)+1) = -1;
                    T(celli(jj), celli(jj)-1) = -1;
            endswitch 
        endfor
    endfor

    D = struct();

endfunction

%!test
%!  spec  = [0 1 1 1 0; 2 0 0 0 0; 0 2 0 0 0; 0 0 2 0 0; 0 0 0 2 0]
%!  sregm = mksmth(spec);
%!  disp(sregm);
%!  error("boo")
