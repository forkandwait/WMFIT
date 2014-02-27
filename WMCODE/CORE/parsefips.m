function [st, cnty] = parsefips(infips)
	## parses fips integer into state and county
		 
    cnty = mod(infips, 1000);
    st = (infips - cnty) ./ 1000;
endfunction
