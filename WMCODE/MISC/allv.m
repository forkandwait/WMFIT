function S = allv()
    ## Collect all variables of the caller into a struct, usually for
    ## debugging.

    S = struct();
    for vn = evalin("caller", "who()")'
        nm = vn{1};
        S = setfield(S, nm, evalin("caller", nm));
    endfor
endfunction
