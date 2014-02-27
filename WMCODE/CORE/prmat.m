function out = prmat(A, doprint=1)
	## print transition matrix

	out = "";

	## get skeleton and schematic stuff
	matstuff = mkmat();

	## format str collector
	fmt = "";
	cfmts = { '%8i' '%8.4f'};

	## make age labels
	out = [out sprintf(repmat("%8.2i", 36), [0:5:85 0:5:85]) "\n"];

	## format each row with specialized format
	for rowi = 1:36
		ni = (matstuff.tmat(rowi,:) >= 1) + 1;
		fmt = [strcat(cfmts(ni){:}) "\n"];
		out = [out sprintf(fmt, A(rowi, :))];
	endfor

	if doprint == 1
		printf(out);
	end

endfunction
