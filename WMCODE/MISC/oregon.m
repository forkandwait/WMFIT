function out =  oregon (P, st="41", isinter=0, varmult=1)
    ## look at oregon

    ## defaults for bands procedure
    _simn = 48;
    _fwdn = 8;
    _varmult = 1.0;

    out = cell();
    fooi = 1;
    for cntyi = (strmatch (st, P.r_geo))'
        _fcstname = P.r_geo{cntyi};

        inp = P.pop{cntyi};
        totalp = sum(inp)(end);

        out{fooi} = bands(inp, _simn, _fwdn, varmult, isinter, _fcstname);

        if isinter
            input(sprintf("I love Oregon! (%s). %i.", _fcstname, totalp));  
        else
            disp(sprintf("I love Oregon! (%s). %i.", _fcstname, totalp));
        endif
        fooi++;
    endfor
endfunction
