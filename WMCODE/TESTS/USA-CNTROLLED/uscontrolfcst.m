function [fstr] = uscontrolfcst (alldata, fcstn=2, strtn=3, stopn=9, colleger=1.4)
		 
	## usage: [fstr] = uscontrolfcst (inpop, fcstn=2, strtn=3, stopn=9, colleger=1.4))
	## 
	## Returns a structure with a whole lot of forecast data

	fstr = struct();
	fstr.input.fcstn = fcstn;
	fstr.input.strtn = strtn;
	fstr.input.stopn = stopn;
	fstr.input.colleger = colleger;

	## get base data at USA, ST, CNTY
	[fstr.input.usa fstr.input.st fstr.input.cnty fstr.input.stfips fstr.input.cntyfips fstr.input.cntystinfo ] = ustotals(alldata);
	usa  = fstr.input.usa(:, strtn:stopn);
	st   = fstr.input.st(:, strtn:stopn, :);
	cnty = fstr.input.cnty(:, strtn:stopn, :);

	## get transition matrix for the US and do a forecast
	[fstr.usa.A, fstr.usa.obj, fstr.usa.info, fstr.usa.lambdas] = wmh(usa);
	fstr.usa.fcst = nan(36, fcstn);
	for yri = 1:fcstn
		fstr.usa.fcst(:,yri) = fstr.usa.A^yri * usa(:,end);
	endfor

	## get transition matrices and NON-controlled forecasts for each state
	staten = (size(st)(3));
	for sti = 1:staten
		[fstr.state.A(:,:,sti) fstr.state.obj(sti) fstr.state.info(sti) fstr.state.lambdas{sti} ] = wmh(st(:,:,sti));
		for yri = 1:fcstn
			fstr.state.fcst_nc(:,yri,sti) = fstr.state.A(:,:,sti)^yri * st(:,end,sti);
		endfor
	endfor

	## control states to usa with wm_alloc
	fstr.state.fcst = round(wm_alloc(fstr.usa.fcst, fstr.state.fcst_nc));

	## get county NON-controlledx forecasts
	cntyn = (size(cnty)(3));
	for cntyi = 1:cntyn
		_cohortr = cnty(2:18,end,cntyi)./cnty(1:17,end,cntyi);
		if _cohortr(4) >= colleger
			warning ("college county! cntyi=%i", cntyi) ;
			_cohortr
			fstr.cnty.modeltype{cntyi} = "college";
		else 
			fstr.cnty.modeltype{cntyi} = "generic";
		endif
		#[fstr.cnty.A(:,:,cntyi) fstr.cnty.obj(cntyi) fstr.cnty.info(cntyi) fstr.cnty.lambdas(:,cntyi)] = wmh(cnty(:,:,cntyi), [], fstr.cnty.modeltype{cntyi});
		[fstr.cnty.A(:,:,cntyi)  fstr.cnty.obj(cntyi) fstr.cnty.info(cntyi) fstr.cnty.lambdas{cntyi}] = wmh(cnty(:,:,cntyi), [], fstr.cnty.modeltype{cntyi});
		for yri = 1:fcstn
			_fcstx = fstr.cnty.A(:,:,cntyi)^yri * cnty(:,end,cntyi);
			_fcstx(_fcstx <= 0) = 0; # keep it from going neg
			fstr.cnty.fcst_nc(:, yri, cntyi) = _fcstx; 
		endfor
	endfor

	## for each state, control its county's to it.
	[_ _ _stind] = unique(fstr.input.stfips(1:36:end));  # determine county indices
	for sti = 1:staten 
		fstr.input.instatei = find(_stind == sti);
		currcntydata = fstr.cnty.fcst_nc(:,:,fstr.input.instatei);
		fstr.cnty.fcst(:, :, fstr.input.instatei) = round(wm_alloc(fstr.state.fcst(:,:,sti), currcntydata)); 
	endfor
endfunction
