function [errfoo out] = do_gma(P, G, isinter=0)
		 
	for cntyi = 1:39
		cntyi
		fips = sprintf("530%2.2i", cntyi*2-1)
		if any (cntyi == [19 37 38])
			constr = mkconwm("college");
            warning ("using college constraints. %i.", cntyi);
		else 
			constr = mkconwm("generic");
		endif
		out(cntyi) = singlecnty(fips, constr, P, G.pop(:,cntyi));
		printf("MAPE me  : %6.6f, them: %6.6f\n", out(cntyi).ERR.MAPE,   out(cntyi).CERR.MAPE);
		printf("MedAPE me: %6.6f, them: %6.6f\n", out(cntyi).ERR.MedAPE, out(cntyi).CERR.MedAPE);
		printf("KLD me   : %6.6f, them: %6.6f\n", out(cntyi).ERR.KLD,    out(cntyi).CERR.KLD);
		
		errfoo.perrwm(:,cntyi)   = abs(out(cntyi).ERR.per);
		errfoo.perrcomp(:,cntyi) = abs(out(cntyi).CERR.per);
		errfoo.kldwm(cntyi)      = out(cntyi).ERR.KLD;
		errfoo.kldcomp(cntyi)    = out(cntyi).CERR.KLD;

        if isinter
           	input("sdlkfj");
        endif

	endfor

    disp([quantile(errfoo.perrwm(:), [0 .15 .5 .85 1]), quantile(errfoo.perrcomp(:), [0 .15 .5 .85 1])]);    
    disp([quantile(errfoo.kldwm(:), [0 .15 .5 .85 1]), quantile(errfoo.kldcomp(:), [0 .15 .5 .85 1])]);
    
endfunction
