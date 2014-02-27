function prcurr(H, subdir='')
    ## print figure based on itnernal name, converting weird characters to  
    fname = strcat(subdir, regexprep(get(H, 'name'), "[\\W]", "_"), ".png");
    print(H, fname);
endfunction
