function out = allcounties(P)
    ## Do all counties with defaults, separate by state

    for _st = States()'
        disp(_st{3}); 
        st_data = oregon(P, _st{1});
        prstate(st_data, _st{2});
        
    endfor
endfunction
