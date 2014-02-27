function out = Pop ()
    ## Returns a struct with USA county data.
         
    _dir = fileparts(mfilename("fullpathext"));
    x = dlmread(fullfile(_dir, "../../DATA/USCOUNTIES/allcntysdata.txt"), "\t");

    out.pop = mat2cell(x(:,4:end), repmat(36, rows(x)/36, 1), 9); 

    ## fill out fips codes
    out.r_geo = cell(rows(x)/36, 1);
    for ii = 1:length(out.r_geo)
        fi = ii*36-35;
        out.r_geo{ii,1} = sprintf("%5.5i", x(fi,1));
    endfor
    

    ## calculate state totals



    out.description = {"All county age and sex, vintage Jan 1, 2102 (original dissertation working data). 1970 to 2010, 5 year.  Caveats apply."};

endfunction
