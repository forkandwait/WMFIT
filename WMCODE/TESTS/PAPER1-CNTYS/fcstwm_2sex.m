function [out] = fcstwm_2sex(datapop, startin, endin, endf, label, strtyr, do_disp=0, do_save=1, colleger=1.4)
    ## [out] = fcstwm_2sex(datapop, startin=3, endin=3, endf=6, label, do_disp, do_save)

    ## NOTE: CURRENTLY 2011-10-19 DOING A TEN YEAR BACK TO NOW FORECAST FOR DISPLAY AND PAPER WRITING.  HACK
    endemp = columns(datapop);

	## choose college classification (only college or generic so far)
	cohratio = reshape(datapop([2:18 20:36],end) ./ datapop([1:17 19:35],(end-1)), [], 2);
	if any(cohratio(4,:) >= colleger) ## && cohratio(5,:) <= 1/colleger)
		modeltype = "college";
		constr = mkconwm(modeltype);
	else
		modeltype="generic";
		constr = mkconwm(modeltype); 
	end
	if ~(strcmp(modeltype, "generic"))
		warning("%s, %s:  using \"%s\" model!", label, mfilename, modeltype);
	endif

    ## save out input values
    out.datapop = datapop;

    ## calculate year values based on step values
    endempyr = strtyr + 5*(endemp-startin);
    endinputyr = strtyr + 5*(endin-1);
    endfcstyr = endempyr + 5*endf;
    
    ## generate the transition matrix, and a testing version with two steps lopped off
    [out.A, _stuff] = wmh(datapop(:,startin:endin), constr);
    out.obj = _stuff.obj;
    out.info = _stuff.info;
    out.lambda = _stuff.lambda;

    ## generate extrapolation matrices and test data
    out.Atest10 = wmh(datapop(:,startin:(endin-2)), constr);
    out.f_test10 = out.Atest10^2 * datapop(:,endin-2);
	if any(isnan([out.Atest10(:); out.f_test10(:)]))
	   warning("%s: nan in test forecast. %s.", mfilename, label);
	endif

	## generate H-P matrix from last pops (hopefully census) and forecasts out to limit
	out.hptrm = hamp(datapop(:,(end-3)), datapop(:,(end-2), 1));
    for fi = 1:4
		out.hpfcst(:,fi) = ((out.hptrm)^fi)*datapop(:, (end-2));
    end

    ## forecast out endf steps from "endin" (might overlap with
    ## empirical input). "Afcst" for matrix "A"
    for fi = 1:endf
	  out.Afcst(:,fi) = (out.A)^fi*datapop(:, endin);
    end

    ## calculate linear fit to totals and fcst
    out.total = sum(datapop(:,startin:endin));
    [p s mu] = polyfit(startin:endin, out.total, 1);
    [out.tfcst dy] = polyval(p, startin:(endemp+endf), s, mu);
    out.tfcst = [out.tfcst'-dy', out.tfcst', out.tfcst'+dy'];

    if do_disp
	 
	  ############################################################
	  ## Make pictures!
	  ############################################################
	  
	  ## various interpolated forecasts for check
	  f_12 = out.A*datapop(:,startin);     # datapop(:,1) to one step ahead
	  f_2end = out.A^2*datapop(:,endin-2); # datapop(:,end-1) to end

	  ## HP and Woods and empirical 
	  figure(4);
	  set(gca, 'colororder', [ .80 .80 .80; 0 0 0; 0 0 1; 1.0 0 0 ]);
	  plot(0:5:85, [datapop(1:18, (end-2)), datapop(1:18, (end)), out.hpfcst(1:18,2), out.f_test10(1:18)]);
	  legend({'2000 p', '2010 p', '2010 hp', '2010 a'});

	  ## males
	  figure(1); 

	  ## forecast one step ahead from first datapop, comparison
	  s1 = subplot(3,1,1);
	  set(gca, 'colororder', [ 1 0 0; 0 0 0; .85 .85 .85 ]);
	  plot(0:5:85, [f_12(1:18), datapop(1:18,startin+1), datapop(1:18,startin)]);
	  legend(sprintf("%i K (model 1stp)", strtyr+5),
		   sprintf("%i K (emp)", strtyr+5), 
		   sprintf("%i K (emp)", strtyr),
		   'location', 'northeast'
		  );  
	  _y1 = ylim();

	  ## forecast to end of data from jumpoff INTERPOLATED, comparison
	  s2 = subplot(3,1,2);
	  set(gca, 'colororder', [ 1 0 0; 0 0 0; .85 .85 .85; ]);
	  plot(0:5:85, [f_2end(1:18), datapop(1:18,endin), datapop(1:18,startin)]);
	  legend(sprintf("%i K (model 2stp interp)", endempyr), 
		   sprintf("%i K (emp)", endempyr),
		   sprintf("%i K (emp)", strtyr)
		   ,'location', 'northeast'
		  );  
	  _y2 = ylim();

	  ## forecast to end of data from jumpoff EXTRAPOLATED, comparison
	  s3 = subplot(3,1,3);
	  set(gca, 'colororder', [ 1 0 0; 0 0 0; .85 .85 .85; ]);
	  plot(0:5:85, [out.f_test10(1:18), datapop(1:18,endin), datapop(1:18,endin-2)]);
	  legend(sprintf("%i K (model 2stp extrap)", endempyr), 
		   sprintf("%i K (emp)", endempyr),
		   sprintf("%i K (emp)", endempyr-10)
		   ,'location', 'northeast'
		  );  
	  _y3 = ylim();

	  #{
	  ## forecast to end of data from jumpoff EXTRAPOLATED, comparison
	  s4 = subplot(4,1,4);
	  set(gca, 'colororder', [ 1 0 0; 0 0 0; .85 .85 .85; ]);
	  plot(0:5:85, [(out.f_test20(1:18)), (datapop(1:18,endin)), (datapop(1:18,endin-4))]);
	  legend(sprintf("%i K (model 4stp extrap)", endempyr), 
		   sprintf("%i K (emp)", endempyr),
		   sprintf("%i K (emp)", endempyr-20)
		   ,'location', 'northeast'
		  );  
	  set (gca, 'xcolor', [.8 .8 .8], 'ycolor', [.8 .8 .8]);
	  grid on;
	  _y4 = ylim();
	  #}

	  ## ylimits
	  _y = max([_y3; _y2; _y1]);
	  set(s1, 'ylim', _y);
	  set(s2, 'ylim', _y);
	  set(s2, 'ylim', _y);
	  	  
  	  ## females
	  figure(2);

	  ## forecast one step ahead from first datapop, comparison
	  s1 = subplot(3,1,1);
	  set(gca, 'colororder', [ 1 0 0; 0 0 0; .85 .85 .85 ]);
	  plot(0:5:85, [f_12(19:36), datapop(19:36, startin+1), datapop(19:36, startin)]);
	  legend(sprintf("%i K (model 1stp)", strtyr+5),
		   sprintf("%i K (emp)", strtyr+5), 
		   sprintf("%i K (emp)", strtyr)
		   ,'location', 'northeast'
		  );  
	  _y1 = ylim();

	  ## forecast to end of data from jumpoff, comparison
	  s2 = subplot(3,1,2);
	  set(gca, 'colororder', [ 1 0 0; 0 0 0; .85 .85 .85; ]);
	  plot(0:5:85, [f_2end(19:36), datapop(19:36,endin), datapop(19:36, startin)]);
	  legend(sprintf("%i K (model 2stp interp)", endempyr), 
		   sprintf("%i K (emp)", endempyr),
		   sprintf("%i K (emp)", strtyr)
		   ,'location', 'northeast'
		  );  
	  _y2 = ylim();

	  ## forecast to end of data from jumpoff EXTRAPOLATED, comparison
	  s3 = subplot(3,1,3);
	  set(gca, 'colororder', [ 1 0 0; 0 0 0; .85 .85 .85;  ]);
	  plot(0:5:85, [out.f_test10(19:end), datapop(19:end,endin), datapop(19:end,endin-2)]);
	  legend(sprintf("%i K (model 2stp extrap)", endempyr), 
		   sprintf("%i K (emp)", endempyr),
		   sprintf("%i K (emp)", endempyr-10)
		   ,'location', 'northeast'
		  );  
	  _y3 = ylim();

	  #{
	  ## forecast to end of data from jumpoff EXTRAPOLATED, comparison
	  s4 = subplot(4,1,4);
	  set(gca, 'colororder', [ 1 0 0; 0 0 0; .85 .85 .85;  ]);
	  plot(0:5:85, [out.f_test20(19:end), datapop(19:end,endin), datapop(19:end,endin-4)]);
	  legend(sprintf("%i K (model 4stp extrap)", endempyr), 
		   sprintf("%i K (emp)", endempyr),
		   sprintf("%i K (emp)", endempyr-20)
		   ,'location', 'northeast'
		  );  
	  grid on;
	  set (gca, 'xcolor', [.8 .8 .8], 'ycolor', [.8 .8 .8]);
	  _y4 = ylim();
	  #}

	  ## ylimits
	  _y = max([_y3; _y2; _y1]);
	  set(s1, 'ylim', [0 _y(2)]);
	  set(s2, 'ylim', [0 _y(2)]);
	  set(s3, 'ylim', [0 _y(2)]);
	  

	  ############################################################
	  ## project endf steps forward, mf separately
	  ############################################################
	  figure(3);

	  ##
	  s1 = subplot(3,1,1);
	  set(gca, 'colororder', [.8 .8 .8; 1 0 0]);
	  plot(0:5:85, [datapop(1:18,end), out.Afcst(1:18,end)]); 
	  legend(sprintf("%i (emp, M)", endempyr),
		   sprintf("%i (model, M, jump=%i)", endfcstyr, endempyr)
		  ); 
	  _y1 = ylim();
	  

	  s2 = subplot(3,1,2);
	  set(gca, 'colororder', [.8 .8 .8; 1 0 0]);
	  plot(0:5:85, [datapop(19:end,end), out.Afcst(19:end,end)]); 
	  legend(sprintf("%i (emp, F)", endempyr),
		   sprintf("%i (model, F, jump=%i)", endfcstyr, endempyr)
		  ); 
	  _y2 = ylim();
	  
	  _y = max([_y2; _y1]);
	  set(s1, 'ylim', [0 _y(2)]);
	  set(s2, 'ylim', [0 _y(2)]);

	  ## totals
	  s3 = subplot(3,1,3);	  
	  set(gca, 'colororder', [.87 .87 .87; .87 .87 .87; .87 .87 .87]);
	  hold off;
	  plot(strtyr:5:endfcstyr, out.tfcst);
	  hold on;
	  plot(strtyr:5:endempyr, out.total, 'k');
	  plot((endempyr+5):5:endfcstyr, sum(out.Afcst), 'r');
	  hold off;

	  if 0
		  for _pi = 1:4
			  refresh(_pi);
		  endfor
	  endif

	  ############################################################
	  ## print matrix nicely   
	  ############################################################
	  printf("\n\n                                    %s\n", label);
	  printf("info = %i, solveiter = %i\n", out.info.info, out.info.solveiter);
	  printf("\tstart Year: %i, \n\tend input year: %i, \n\tend training year: %i \n\tend fcst year: %i\n", 
		   strtyr, endempyr, endinputyr, endfcstyr);
	  printf("\tN steps input: %i, N steps data: %i, N steps past jump: %i\n", endin, endemp, endf);
	  tfr = sum(out.A([1 19],21:29)(:)); 
	  printf("Error pctiles (10 yr): [%2.2f%% %2.2f%% %2.2f%% %2.2f%% %2.2f%%] (0, .25, .50, .75, 1)\n", 
			 d_err(datapop(:,end), out.f_test10));
	  printf("Total \"fertility\": %f\n", tfr);
	  printf("Total pop at jumpoff: %i\n", sum(datapop(:,end)));
	  ## print the matrix (tricky, a little)
	  printf("  ");
	  printf(repmat("      %2.2i", 36), [0:5:85 0:5:85]);
	  printf("\n"); 
	  _fmt = repmat("%7.4f ", 1, columns(out.A));
	  printf(["%1.2i " _fmt "\n"], [[0:5:85 0:5:85]', out.A]');
    end;

end;

function d_err =  d_err(emp, pred)
	## calculates an error mesure for an emp and a pred vector.  Currently MedianAPE
	perrv = round(100*abs(pred(:)-emp(:))./emp(:));
	d_err = quantile(perrv, [0 .25 .50 .75 1.00], 1, 3);
end;
