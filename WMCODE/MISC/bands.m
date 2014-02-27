function out = bands(pop, fcstname=[], varmult=1.0,  simn=124, fwdn=2)
    ## generate forecast bands using residual-perturbed input

    ## real forecast
    _p = pop(:,end-2:end);
    out.fcst = runpop(_p, simn, fwdn, [.05 .95], varmult);
    out.fcst.fcstname = fcstname;    
    
    ## xvalid forecast
    _p = pop(:,end-4:end-2); 
    out.xfcst = runpop(_p, simn, fwdn, [.05 .95], varmult); 
    out.xfcst.fcstname = fcstname;
    
    ## plot both
    runplot(out.fcst, [], 1, sprintf("fcst-%s", fcstname));
    runplot(out.xfcst, out.fcst.pop(:,end), 3, sprintf("xvalid-%s", fcstname)); 

endfunction

function out = runpop (pop, simn, fwdn, ptile=[.05 .95], varmult=1.0)

    ## variance fudge factor
    if varmult ~= 1.0 
       ;
    endif

    ## start with Leslie matrix 
    A = wmh(pop);

    ## create modeled pop
    mpop = nan(size(pop));
    mpop(:,1) = pop(:,1);
    for yri = 2:columns(pop)
        mpop(:,yri) = A^(yri-1) * mpop(:,1);
    endfor

    ## calculate residuals, mean, cov. Use eig trick to force positive
    ## semidef (XXX) 
    rpop = pop - mpop;
    rpopt = rpop(:,2:end)';
    rpopm = mean(rpopt);
    [u s] = eig(cov(rpopt));
    rpops = u*abs(s)*u';

    ## simulate lots of populations and Leslie mats
    simP = nan([size(pop) simn]);
    simA = nan([36 36 simn]);
    for simi = 1:simn
        simi
        ## generate simulated residuals
        rpopsim(:,:,simi) = (mvnrnd(rpopm, varmult*rpops, columns(pop)-1))';

        ## add residuals to pop
        simP(:,:,simi) = pop + [zeros(36,1), rpopsim(:,:,simi)];
        if any(simP(:)) < 0
            error("%s negative pop in simulation.", mfilename);
        endif

        ## derive WM Leslie matrix from simulated pop
        simA(:,:,simi) = wmh(simP(:,:,simi));
    endfor

    ## forecast each population fwdn steps
    simFcst = nan(36, fwdn, simn);
    for simi = 1:simn
        for stepi = 1:fwdn
            simFcst(:, stepi, simi) = simA(:,:,simi)^stepi * pop(:,end);
        endfor
    endfor

    ## summarize with quantiles -- THESE GIVE UPPER AND LOWER
    simFcstq = squeeze(quantile(simFcst, ptile, 3)); 

    ## Do a forecast with real A for comparison
    Fcst = nan(36, fwdn);
    for stepi = 1:fwdn
        Fcst(:,stepi) = A^stepi * pop(:,end);
    endfor 

    ## output all variables
    out = allv();
endfunction


function out = runplot(fcst, xpop=[], strtfig=1, fcstname=[])

    ## colors and stuff for plots
    tracecol = [.9 .9 .9];
    ptilecol = [.75 .75 .75];
    isvis = ifelse(isempty(fcstname), true, false);

    ## residuals histogram
    fh(1) = figure(strtfig, 'Name', 'Residual Histogram', 'Visible', isvis);
    foo =fcst.rpop(:);
    hist(foo, 48);

    ## forecasts 
    fh(2) = figure(strtfig+1, 
                   'Name', ifelse(isempty(xpop), 'Forecasts', 'Xvalid Forecasts')); 

    ## set legend for 
    _plegstr = ifelse(isempty(xpop), "Jump-off Pop", "X-valid End Pop");
    
    ## males
    _s(1) = subplot(2,1,1); 
    _h(1) = plot(0:5:85, fcst.Fcst(1:18,end), '-^k', 'linewidth', 2.0); # forecast from input
    hold on;
    if isempty(xpop)
        _h(2) = plot(0:5:85, fcst.pop(1:18,end), '-ok', 'linewidth', 2.0); # end pop 
    else 
        _h(2) = plot(0:5:85, xpop(1:18,end), '-ok', 'linewidth', 2.0); # end pop for x valid 
    endif
    _h([3 4]) = plot(0:5:85, squeeze(fcst.simFcstq(1:18,end,:)), ':k', 'linewidth', 2.0);
    _ylim(1,:) = ylim();
    legend(_h(1:2), "Forecasted Population", _plegstr, sprintf("%i%%", 100*fcst.ptile(1)), sprintf("%i%%", 100*fcst.ptile(2)));
    set(gca, 'xtick', 0:10:90);
    grid on;
    hold off;

    ## females
    _s(2) = subplot(2,1,2);
    hold on;
    _h(1) = plot(0:5:85, fcst.Fcst(19:36,end), '-^k', 'linewidth', 2.0);
    if isempty(xpop) 
        _h(2) = plot(0:5:85, fcst.pop(19:36,end), '-ok', 'linewidth', 2.0); # end pop
    else 
        _h(2) = plot(0:5:85, xpop(19:36,end), '-ok', 'linewidth', 2.0); # end pop for x valid
    endif 
    _h([3 4]) = plot(0:5:85, squeeze(fcst.simFcstq(19:36,end,:)), ':k', 'linewidth', 2.0);
    _ylim(2,:) = ylim(); 
    legend(_h(1:2), "Forecasted Population", _plegstr, sprintf("%i%%", 100*fcst.ptile(1)), sprintf("%i%%", 100*fcst.ptile(2)));
    set(gca, 'xtick', 0:10:90);
    grid on;
    hold off;

    ylim(_s(1), [0 max(_ylim(:,2))]);
    ylim(_s(2), [0 max(_ylim(:,2))]);

    if ~isempty(fcstname)
       warning("bands: printing. Dir: %s.", pwd);
       print(fh(1), sprintf("%s-histr.png", fcstname));
       print(fh(2), sprintf("%s-fcst.png", fcstname));
    else 
        warning("bands: not printing");
    endif

    close (fh(2));

endfunction;

function out = multipops(inpop);
    ## calculate empirical bands by doing successive 10 year step
    ## aheads XXX doesn't work right yet...

    error("function doesn't work");
    stopyr = (columns(inpop)-4);
    for yi = 1:stopyr
        A = wmh(inpop(:,yi:yi+2));
        xfcst(:,yi) = A^2*inpop(:,yi+2);
        out.xfcstr(:,yi) = inpop(:, yi+4) - xfcst(:,yi);
    endfor
    out.pe = out.xfcstr ./ inpop(:,5:end);
    round(median(out.pe,2)*100);
    figure(5); hist(out.xfcstr(:), 32);
    
endfunction
