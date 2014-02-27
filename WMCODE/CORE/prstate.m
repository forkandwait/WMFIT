function out = prstate(st_bands, fname)
    ## prints state data
    ## prstate('01', counties, "AL.txt") 
		 
    ## real forecast
	[fid fmsg] = fopen(sprintf("%s.txt", fname), 'w');
	if fid <= 0
		error(sprintf("%s: %s", mfilename, fmsg));
	endif
    loopcountys(fid, st_bands, 0);
    fclose(fid);

    ## xvalid forecast
	[fid fmsg] = fopen(sprintf("%s-xvalid.txt", fname), 'w');
	if fid <= 0
		error(sprintf("%s: %s", mfilename, fmsg));
	endif
    loopcountys(fid, st_bands, 1);
    fclose(fid);

endfunction

function out = loopcountys (fid, st_bands, isxvalid=0)
         
    fprintf(fid, "Copyright %s, W. Webb Sprague.\n\n", datestr(date,10));
    fprintf(fid, "By using these forecasts, you agree to cite http://arxiv.org/abs/1203.2313 in all direct or derived works.\n\n\n");
	for cntyi = 1:columns(st_bands)
        cnty = st_bands{1, cntyi};
        if isempty(cnty) || isempty(cnty.fcst) || isempty(cnty.xfcst)
           continue
        endif 
        fcst = ifelse(isxvalid, cnty.fcst, cnty.xfcst);
        
        ## county heading
        fprintf(fid, "\n***************************************\n");
        fprintf(fid, "    County %s\n\n", fcst.fcstname);

        ## low forecast
        fprintf(fid, "Low Forecast\n");
        fdisp(fid, round(fcst.simFcstq(:,:,1)));
        fprintf(fid, "\n\n");

        ## med forecast
        fprintf(fid, "Med Forecast\n");
        fdisp(fid, round(fcst.Fcst));
        fprintf(fid, "\n\n");

        ## high forecast
        fprintf(fid, "High Forecast\n");
        fdisp(fid, round(fcst.simFcstq(:,:,2)));
        fprintf(fid, "\n\n");

        ## input data
        fprintf(fid, "Input Data\n");
        fdisp(fid, round(fcst.pop));
        fprintf(fid, "\n\n");
    endfor
    out = 1;
endfunction

%!demo 
%!  prstate('01', counties, "foo.txt")
