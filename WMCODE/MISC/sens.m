function [senmat elmat wmstuff] = sens(A)
    ## calculate sensitivities and elasticities for A, from Caswell
    ## chapter 9. 

    warning("%s using pinv() because of ill-conditioned A matrices.", mfilename());
         
    [W, d] = eig(A);
    d = diag(d);
    imax = find(d==max(d));
    V = conj(pinv(W));
    w = W(:, imax);
    v = real(V(imax,:))';

    senmat = v* w';
    elmat = (senmat .* A) / max(d);

    ## parse out if 36 x 36 WM 
    if all(size(A) == [36 36])
       _m = mkmat();
       wmstuff.mig  = [elmat(_m.tmat==1), elmat(_m.tmat==4)];
       wmstuff.surv = [elmat(_m.tmat==2), elmat(_m.tmat==5)];
       wmstuff.surv = [wmstuff.surv; 0 0];
       wmstuff.fert = [elmat(_m.tmat==3), elmat(_m.tmat==6)];
       wmstuff.fert = [0 0 ; wmstuff.fert; zeros(9, 2)];
    endif
endfunction
