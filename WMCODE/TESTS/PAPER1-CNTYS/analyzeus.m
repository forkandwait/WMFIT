function [out err_rank] = analyzeus(fdata, subsetf="none")
    ## analyzeus(fdata, subsetf="none") (error percentiles etc)

    ## subset counties
	out.subsetf = subsetf;
    switch subsetf
		case '<50k'
			MINCNTY	    = 50000;
			subi	    = find(sum(fdata.emp10yr) < MINCNTY);
		case '>=50k'
			MINCNTY	    = 50000;
			subi	    = find(sum(fdata.emp10yr) >= MINCNTY);
		case '>=20k'
			MINCNTY	    = 20000;
			subi	    = find(sum(fdata.emp10yr) >= MINCNTY);
		case 'wa'
			subi        = 2932:2970;
		case "none"
			subi = 1:length(fdata.emp10yr);
		otherwise
			error("%s: unknown subset option: \"%s\"", mfilename, subsetf);
    endswitch

    ##     set intervals and histogram edges for robust analysis
    _q  = [0.50 0.80 0.90 0.95 0.975 0.99 1.00];
    _bw = 0.025;
    _x  = 0:_bw:1;

	## set up useful variables
    emp10yr	    = fdata.emp10yr(:, subi);
    emp00yr	    = fdata.emp00yr(:, subi);
    f10yr	    = fdata.f10yr(:, subi);
	hpf10yr     = fdata.hpf10yr(:, subi);
    infocnty    = fdata.infocnty(subi);
    labels	    = fdata.labels(subi);
	
    ## add a single person to any zero cells, to avoid Infs
    emp10yr(emp10yr == 0) = 1;
    f10yr(f10yr == 0)     = 1;

	## start delete process
	badi = [];

	## delete NaN HP forecasts (XXX ???)
	_badi = find(any(isnan([hpf10yr])));
	badi = [badi _badi];
	if ~isempty(_badi)
		showacky("NAN HP fcsts", _badi, labels);
	end

	## delete NaN W forecasts (XXX ???)
	_badi = find(any(isnan([f10yr])));
	badi = [badi _badi];
	if ~isempty(_badi)
		showacky("NAN W fcsts", _badi, labels);
	end

    ## delete non convergent forecasts
    optinfo = nan(length(infocnty),1);
    for ii = 1:length(infocnty)
	  optinfo(ii) = infocnty(ii).info.info;
    end
    _badi = find(optinfo ~= 0)';
	badi = [badi _badi];
	if ~isempty(_badi)
		showacky("Bad optimizations", _badi, labels);
	end

    ## deal with wacky large results (should be nonconvergent)
    _badi = find(any(abs(f10yr ./ emp10yr) > 100));
	badi = [_badi badi];
	if ~isempty(_badi)
		showacky("Wacky large populations", _badi, labels);
	end

	## delete them all.
	printf("\nDELETING %i WEIRD FORECASTS.\n ", length(badi));
	emp10yr(:,badi) = [];
	emp00yr(:,badi) = [];
	f10yr(:,badi)   = [];
	hpf10yr(:,badi) = [];
	infocnty(badi)  = [];
	labels(badi)    = [];

	## output the deleted versions
	out.emp10yr = emp10yr;
	out.emp00yr = emp00yr;
	out.f10yr = f10yr;
	out.hpf10yr = hpf10yr;
	out.infocnty = infocnty;
	out.labels = labels;


    ## interesting indices
    m   = 1:18; 
    f   = 19:36;

    ## male and female
    emp10yrm	= emp10yr(m, :);
    emp10yrf	= emp10yr(f, :);
    f10yrm	    = f10yr(m, :);
    f10yrf	    = f10yr(f, :);
	hpf10yrm    = hpf10yr(m, :);
	hpf10yrf    = hpf10yr(f, :);

    ## totals 
    emp10yrmt	= sum(emp10yrm);
    emp10yrft	= sum(emp10yrf);
    f10yrmt		= sum(f10yrm);
    f10yrft		= sum(f10yrf);

    ## age structure
    emp10yrms	= bsxfun('rdivide', emp10yrm, emp10yrmt);
    emp10yrfs	= bsxfun('rdivide', emp10yrf, emp10yrft);
    f10yrms	= bsxfun('rdivide', f10yrm, f10yrmt);
    f10yrfs	= bsxfun('rdivide', f10yrf, f10yrft);

    ## age count normed to total
    f10yrmn     = bsxfun('times', emp10yrmt, f10yrms);
    f10yrfn     = bsxfun('times', emp10yrft, f10yrfs);

    ## abs percentage error by agesex and county
    dm		= emp10yrm - f10yrm;
    df		= emp10yrf - f10yrf;
    pferrm	= abs(dm ./ emp10yrm);
    pferrf	= abs(df ./ emp10yrf);

	## compute KL Divergence (on proportion) for all counties 
	_klerr = nan(length(emp10yr),1);
	for ii = 1:length(emp10yr)
		_klerr(ii) = kldist(emp10yr(:,ii), f10yr(:,ii)); 
	end
	out.klerr = _klerr;
	[klerr _klerri] = sort(_klerr(:), "descend");

	## total percentage error for each county, sorted and associated
	## with FIPS codes
	_tpop = round(sum(emp10yr)');
	_terr = sum([pferrm; pferrf])';
	out.terr = _terr;
	[_terrs _terri] = sort(_terr, "descend");

	## make an error report
	err_rank = cell(length(_terri), 4);
	err_rank(:,1) = labels(1,:);
	err_rank(:,2) = num2cell(_tpop);
	err_rank(:,3) = num2cell(_terr);
	err_rank(:,4) = num2cell(_terri);
	err_rank(:,5) = num2cell(_klerr);
	err_rank(:,6) = num2cell(_klerri);

	err_rank(:,5) = num2cell(_terri);
	err_rank(:,6) = num2cell(_klerri);

	_err_rank = err_rank(_terri,:);

	_err_rank = err_rank(_klerri,:);

	##     robust SAPP error over all counties
	qterr = quantile(_terr, _q, 1, 3);
	[_hmm, _hmm2, qterri] = intersect(qterr, _terr);
	qterrl = (labels(qterri));
    printf("\nSAPP error over all counties (median = 50%%)\n");
    printf("PERCTILE:"); printf("%5.0i%% ", _q*100); printf("\n");
    printf("FIPS:    "); printf("%6s ", qterrl{:}); printf("\n");
    printf("----------------------------------------------------------\n         ");
    printf("%5.4f ", qterr);
    printf("\n");

	##     robust KL error over all counties
	qklerr = quantile(_klerr, _q, 1, 3);
	[_hmm, _hmm2, qklerri] = intersect(qklerr, _klerr);
	qklerrl = (labels(qklerri));
    printf("\nKL error over all counties (median = 50%%)\n");
    printf("PERCTILE:"); printf("%5i%% ", _q*100); printf("\n");
    printf("FIPS:    "); printf("%6s ", qklerrl{:}); printf("\n");
    printf("---------------------------------------------------------\n         ");
    printf("%4.4f ", qklerr);
    printf("\n");

	## MAPE -- mean absolute percentage error -- total and by age group
	mape = mean ([pferrm, pferrf](:));
	printf("\nMAPE total: %i%%\n", mape*100);
	mape = mean ([pferrm, pferrf], 2);
	printf("\nMAPE by age\n");
	printf("---------\n");
	printf("%2i %2i%%\n", [(0:5:85)', mape*100]');

    ## quantile errors -- uses funky sas method #3 to get findable indexes




    ##     robust percentage errors over all cells 
    allerr = ([pferrm; pferrf])(:);
    qallerr = quantile(allerr, _q, 1, 3);
    printf("\nAbsolute percent quantile errors, all cells (medianAPE = 50%%)\n");
    printf("Age/ Perc"); printf("%4i%% ", _q*100); printf("\n");
    printf("--------------------------------------------------\n         ");
    printf("%4i%% ", qallerr*100);
    printf("\n");

    ##     picture thereof
    figure(1);
    bar(_x, histc(allerr, _x)/length(allerr));
    title("Percent error (over all cells)");
    ylabel("proportion");
    xlabel(sprintf("percent error/100 (e.g. 0.2 = 20%%), binwidth = %2.2f%%", _bw*100));
    grid on;

    ## errors at each age -- absolute
    s = ['M' 'F'];
    qmaerr	= quantile(pferrm, _q, 2, 3); 
    qfaerr	= quantile(pferrf, _q, 2, 3);
    x = [qmaerr; qfaerr] * 100;
    printf("\nM and F absolute percent quantile errors (medianAPE = 50%%)\n");
    printf("Age/ Perc"); printf("%4i%% ", _q*100); printf("\n");
    printf("--------------------------------------------------\n");
    for ii = 1:36
	  printf("%2i %c     ", mod((ii-1)*5, 90), s(floor(ii/19)+1));
	  printf("%4i%% ", x(ii,:));
	  printf("\n");
    end;  

end

function showacky(announce, badi, fiplabels)
	_labs = horzcat(fiplabels(badi)', num2cell(badi)');
	printf('\n%s: \n', announce);
	printf("fips: %s, %i\n", _labs'{:});
end

function d = kldist(epop, fpop)
	## compute kullblack-liebler divergence
	P = epop(:) ./ sum(epop(:));
	Q = fpop(:) ./ sum(fpop(:));
	_d = P .* log(P ./ Q);
	d = abs(sum(_d));
end
	
