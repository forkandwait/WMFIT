function [smpop_c] = wm_alloc(bpop, smpop)
	## controls smpop (3-d: age, county, time) to bigpop (2-d: age, time)
	##
	## [smpop_c] = wm_alloc(bpop, smpop)

	smpop_c = nan(size(smpop));

	for yri = 1:columns(bpop)
		## compute small pop totals using proportional allocation
		## against big pop

		## total big pop across ages this year 1
		_bpop = sum(bpop(:,yri), 1); 

		## age x cnty this year 36 x 51
		_smpop = squeeze(smpop(:,yri,:)); 

		## each county total and proportion via smpop
		_smpopt = sum(_smpop, 1);	# 1x51
		_smpoptpr = (_smpopt ./ sum(_smpopt)); # 

		## reallocated each county totaol
		_smpopta = _smpoptpr * _bpop;

		## feed and run IPF function
		C = _smpopta;
		R = bpop(:,yri);
		T = _smpop;

		smpop_c(:,yri,:) = ipf(T, R, C);

	endfor

endfunction

%!test
%!    smpop = [10 15 10; 10 20 20; 5 5 10];
%!    bpop  = [30; 45; 22];
%!    wm_alloc(bpop, smpop)

%!test
%!    smpop(:,:,1) = [10 15 10; 10 20 20; 5 5 10];
%!    smpop(:,:,) = [10 15 10; 10 20 20; 5 5 10];
%!    bpop  = [30; 45; 22];
%!    wm_alloc(bpop, smpop)
