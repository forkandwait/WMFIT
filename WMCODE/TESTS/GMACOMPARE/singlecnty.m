function out =  singlecnty (fips, constr, P, comp=[], strtp=3, endp=2, doprint=0)
         
    ## blah

    ## prep stuff
    popid = strmatch(fips, P.r_geo);
    if isempty(popid)
       error("Non existent fips: %s", fips);
    endif

    ## handle multi fips (like states)
    if length(popid) >= 2
       pop3 = cat(3, P.pop{popid});
       out.inpop = sum(pop3, 3);
    else
        out.inpop = P.pop{popid};
    endif

    out.endp = endp; 
    out.strtp = strtp;     

    
    ## do the WM thing
    [out.OPT.A, out.OPT.obj, out.OPT.info, out.OPT.lambda] = wmh(out.inpop(:,strtp:end-endp), constr);
    out.pp = round((out.OPT.A^endp) * out.inpop(:,end-endp));

    ## derive errors WM
    out.ERR.diff = out.pp - out.inpop(:,end);
    out.ERR.per  = out.ERR.diff ./ out.inpop(:,end);

    out.ERR.MAPE = mean(abs(out.ERR.per));
    out.ERR.MALPE = mean(out.ERR.per);
    out.ERR.MedAPE = median(abs(out.ERR.per));
    
    pp_prop = out.pp ./ sum(out.pp);
    emp_prop = out.inpop(:,end) ./ sum(out.inpop(:,end));
    out.ERR.KLD = sum(emp_prop .* (log(emp_prop ./ pp_prop)));

	## derive errors for comparison
	if ~isempty(comp)
		out.comp = comp;

		out.CERR.diff = out.comp - out.inpop(:,end);
		out.CERR.per  = out.CERR.diff ./ out.inpop(:,end);
		
		out.CERR.MAPE = mean(abs(out.CERR.per));
		out.CERR.MALPE = mean(out.CERR.per);
		out.CERR.MedAPE = median(abs(out.CERR.per));
		
		comp_prop = out.comp ./ sum(out.comp);
		out.CERR.KLD = sum(emp_prop .* (log(emp_prop ./ comp_prop)));

		## plot against comparison
        ## ... actual numbers ...
		figure(1)
        subplot(2,1,1);
		set(gca, 'colororder', [.90 .90 1.0; .90 .90 .90; 0 0 0; 1 0 0.2; 0 .5 1]);
		plot(0:5:85, [out.inpop(1:18,[end-endp-2, end-endp, end]), out.pp(1:18), out.comp(1:18)], 'o-');
        xlabel ("Male age");
		legend ('jmp-2', 'jmp', 'end', 'wm-fcst', 'comp-fcst');
        subplot(2,1,2);
		set(gca, 'colororder', [.90 .90 1.0; .90 .90 .90; 0 0 0; 1 0 0.2; 0 .5 1]);
		plot(0:5:85, [out.inpop(19:36,[end-endp-2, end-endp, end]), out.pp(19:36), out.comp(19:36)], 'o-');
        xlabel ("Female age");
		legend ('jmp-2', 'jmp', 'end', 'wm-fcst', 'comp-fcst');

        ## ... percentage error
		figure(2)
        subplot(2,1,1);
		set(gca, 'colororder', [1 0 0.2; 0 .5 1]);
		plot(0:5:85, abs([out.ERR.per(1:18), out.CERR.per(1:18)]), 'o-');
        subplot(2,1,2);
		set(gca, 'colororder', [1 0 0.2; 0 .5 1]);
		plot(0:5:85, abs([out.ERR.per(19:36), out.CERR.per(19:36)]), 'o-');

        if doprint 
           print(1, sprintf("%s_m.jpg", fips), '-landscape');
           print(2, sprintf("%s_f.jpg", fips), '-landscape');
        endif
	endif

endfunction;


%!demo
%!    P = Pop;
%!    constr = constr = mkconwm("generic");
%!    out = singlecnty('53001', constr, P);

%!demo
%!    P = Pop;
%!    G = gma2002;
%!    constr = mkconwm("generic");
%!    out = singlecnty("53001", constr, P, G.pop(:,1));

%!demo
%!    P = Pop;
%!    G = gma2002;
%!    constr = mkconwm("generic");
%!    out = singlecnty("53", constr, P, sum(G.pop,2));

