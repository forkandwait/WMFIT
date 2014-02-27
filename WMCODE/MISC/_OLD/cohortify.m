function [outm cyrs] = cohortify(mat, strtyr=[], incyr=5)
    ## create a matrix with input diagonals as output columns
    ## (cohorts), padded with nans.
    ##
    ## optionally returns a matrix of the cohort years calculated by
    ## strtyr representing the first birth cohort avail (1970 with my
    ## data)

    [m n] = size(mat);
    cohn = m + n -1;

    outm = nan(m,cohn);

    for cohi = 1:cohn
        ii = -cohn + n - 1 + cohi;
        cdata = diag(mat,ii);
        [s,~] = size(cdata);
        

        if ii < 0
            rstrt = m - cohi + 1;
            rstp  = rstrt + s - 1; 
            outm(rstrt:rstp,cohi) = cdata;
        elseif ii == 0
            outm(1:s,cohi) = cdata;
        elseif ii > 0
            outm(1:s,cohi) = cdata;
        else
            error("wtf?");
        endif
    endfor

    if ~isempty(strtyr)
       ##keyboard
       cstrtyr = strtyr - incyr*(m-1);
       cstpyr  = strtyr + incyr*(n-1);
       cyrs = (cstrtyr:incyr:cstpyr)';
       cyrs = [cyrs (1:rows(cyrs))'];
    endif
endfunction;

%!demo
%!    foo = cohortify(reshape(1:24,8,[]));

#{
function f(ii, wapop)
    figure(ii);
    mf2t = @(x) x(1:18,:) + x(19:36,:); 
    x = mf2t(wapop (:,:,ii));
    [xc yrs] = cohortify(x, 1970);
    plot(0:5:85, xc(:,10:2:24), '-o');
    legend(num2str(yrs(10:2:24,1)));
    grid on;
endfunction

#}
