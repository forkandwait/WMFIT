function out = rpt_allcntys (fdata, trmats, cpop, colleger=1.4)
    ## make dataset for arxiv
         
    ## get filenames and fids
    _dnow = datestr(now, 30);
    disp(_dnow);
    outf_f = sprintf('allcntysrpt-%s-fcst.txt', _dnow);
    outf_m00 = sprintf('allcntysrpt-%s-mat00.txt', _dnow);
    outf_m10 = sprintf('allcntysrpt-%s-mat10.txt', _dnow);
    outfid_f   = fopen(outf_f, 'w');
    outfid_m00 = fopen(outf_m00, 'w');
    outfid_m10 = fopen(outf_m10, 'w');

    ## column names for forecast file
    fprintf(outfid_f, sprintf("%6.6s %6.6s %6.6s %6.6s %6.6s %6.6s %6.6s %6.6s %6.6s\n",
                            "E2000", "F2010", "E2010", "F2015", "F2020", "F2025", "F2030", "F2035", "F2040"));

    for cntyi = 1:length(fdata.emp00yr)
        disp (fdata.labels{cntyi});

        ## output 2000-2010 transition matrix
        matstr = prmat(trmats(cntyi).Atest10, 0);
        fdisp(outfid_m00, fdata.labels{cntyi});
        fdisp(outfid_m00, matstr);
        fdisp(outfid_m00, '');

        ## get 2010fwd Leslie matrix, based on 1980 to 2010
        datapop = cpop{cntyi};
        cohratio = reshape(datapop([2:18 20:36],end) ./ datapop([1:17 19:35],(end-1)), [], 2);
	    constr = ifelse(any(cohratio(4,:) >= colleger), mkconwm("college"), mkconwm("generic"));
        wm10 = wmh(datapop(:, 3:end), constr);

        ## output 2010fwd matrix
        matstr = prmat(wm10, 0);
        fdisp(outfid_m10, fdata.labels{cntyi});
        fdisp(outfid_m10, matstr);
        fdisp(outfid_m10, '');

        ## do 2010fwd forecast 
        _fcstpop = nan(36,6);
        for stepi = 1:6
            _fcstpop(:, stepi) = (wm10 ^ stepi) * fdata.emp10yr(:, cntyi);
        endfor

        ## output 2010fwd forecast with empirical and test years
        _testpop = horzcat(fdata.emp00yr(:,cntyi), fdata.f10yr(:,cntyi), fdata.emp10yr(:,cntyi));
        fcst = horzcat(_testpop, _fcstpop);
        fdisp(outfid_f, fdata.labels{cntyi});
        fdisp(outfid_f, round(fcst));
        fdisp(outfid_f, '');

    endfor

    ## close file ids
    fclose(outfid_f  );
    fclose(outfid_m00);
    fclose(outfid_m10);

endfunction;
