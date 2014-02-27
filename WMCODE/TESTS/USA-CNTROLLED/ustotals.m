function [us, st, cnty, stfips, cntyfips, cntystinfo] = ustotals(alldata);
	## Take big normalized data and deliver total, state, and cnty
	## arrays
	## 
	## "squeeze (st(:,:,3))" should yield a state matrix (for feeding into wmh)
	## 
	## "us" is ready to roll like it is
	## 

	## get basics 
    [stfips cntyfips] = parsefips(alldata(:,1));
	_alldata = alldata(:, 4:end);
    _yrn = columns(_alldata);

	## get counties 
	[cntyu cntyi cntyj] = unique(alldata(:,1));
	cntyn = length(cntyu);
	cnty = nan(36, _yrn, cntyn);
	for ci = 1:cntyn
		_rng = (1+(ci-1)*36):(ci*36); 
		cnty(:,:,ci) = _alldata(_rng,:);
	endfor

	## Get states by summing respective counties year by year.
	## Trickiness with unique() and accumarray().
    st = nan(36, _yrn, 51);
    [stu sti stj] = unique(stfips);
    for yri = 1:_yrn
        x = accumarray([alldata(:,[3 2]), stj], _alldata(:,yri));
        x = reshape(x(:), 36, []);
        st(:, yri, :) = x;
    endfor
	
	## figure out index relating counties into states
	cntystinfo = [stu(stj(1:36:end)), stj(1:36:end), (1:3120)'];


	## get usa
    us = squeeze(sum(st, 3));

endfunction;
