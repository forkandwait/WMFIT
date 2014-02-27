function [out] = mkmat(mtype="twosex")
	## [out] = mkmat(mtype="twosex")
	##
	## create a transition matrix template with {1,2,3,4,5,6}
    ## for parameter locations, get indices of components
    ## (fert, diag, subd).


    switch mtype
        case "twosex"
            ## prep empty matrix to hold stuff
            _x = zeros(36);

            ## fill pieces with different codes
            _x(sub2ind([36 36], 1:18, 1:18))       = 1; # male diag
            _x(sub2ind([36 36], 2:18, 1:17))       = 2; # male subd
            _x(sub2ind([36 36], ones(1,8), 21:28)) = 3; # male fert
	        
            _x(sub2ind([36 36], 19:36, 19:36))        = 4; # fem diag
            _x(sub2ind([36 36], 20:36, 19:35))        = 5; # fem subd
            _x(sub2ind([36 36], ones(1,8)+18, 21:28)) = 6; # fem fert subd
            
	        ## determine the N by 1 indices of the above places in the matrix
	        nzdiagm    = find(_x == 1);
	        nzsubdm    = find(_x == 2);
	        nzfertm    = find(_x == 3);
	        nzdiagf    = find(_x == 4);
	        nzsubdf    = find(_x == 5);
	        nzfertf    = find(_x == 6);

	        ## make big vector with all indices
	        nz = sort([nzfertm; nzfertf; nzsubdm; nzsubdf; nzdiagm; nzdiagf]);
	        out.nz = nz;
	        
	        ## build indices specific to various components of matrix
	        out._fim = find(ismember (nz, nzfertm));
	        out._fif = find(ismember (nz, nzfertf));
	        out._sim = find(ismember (nz, nzsubdm));
	        out._sif = find(ismember (nz, nzsubdf));
            out._dim = find(ismember (nz, nzdiagm));
            out._dif = find(ismember (nz, nzdiagf)); 

	        ## template matrix for ref
	        out.tmat = _x;
        otherwise
            error("unknown mtype: %s", mtype);
    endswitch
	
endfunction
