function [ fdata, trmats, cpop, _cnty ] = runuscntys_2sex_all(STARTIN=3, ENDIN=9, ENDF=4,
                             savenm="countyfcst", 
                             cntyf="../../../DATA/USCOUNTIES/allcntysdata.txt")
    ## load input data, forecast all counties, save to fdata struct and file.  

    ## XXX runuscntys_2sex_all takes care of makeing our test forecasts magically....

    ## stuff
    _cnty = dlmread (cntyf, '\t');    
    STRTYR  = 1970 + (STARTIN-1)*5;
    printf("starting year: %i\n", STRTYR);
    
    ## make a cell with each county pop (right 3 columns are label stuff)
    cpop = mat2cell(_cnty(:,4:end), repmat(36,  size(_cnty)(1)/36 , 1))';
    
    ## Get county fips 
    _cfips = _cnty(1:36:rows(_cnty),1);
    CFIPS = num2cell(_cfips);
    
    ## do each county
    f10yr = nan(36, length(cpop));
    hpf10yr = nan(36, length(cpop));
    emp10yr = nan(36, length(cpop));
    emp00yr = nan(36, length(cpop));
    labels = cell(1, length(cpop));
    infocnty(length(cpop)).x = 1;
    for ii = 1:length(cpop)
        ## get fips county label
        labels(ii) = sprintf("%5.5i", CFIPS{ii});

        ## ... do a little forecast with full Woods method
        out = fcstwm_2sex(cpop{ii}, STARTIN, ENDIN, ENDF, labels{ii}, STRTYR, 0);
        
        ## save 10 year ahead and empirical
        emp00yr(:,ii)  = out.datapop(:, end-2);
        emp10yr(:,ii) = out.datapop(:, end);
        f10yr(:,ii) = out.f_test10;     

        ## do a H-P forecast to 2010 using 1990-2000 data
        _hptr(:,:,1) = hamp(out.datapop(:,3), out.datapop(:,5), 2);
        _hptr(:,:,2) = hamp(out.datapop(:,5), out.datapop(:,7), 2);
        hptr = mean(_hptr,3);
        hpf10yr(:,ii) = hptr * out.datapop(:,7);

        ## solver info
        infocnty(ii).info = out.info;
        
        printf("%s\n", labels{ii});

        trmats(ii).A = out.A;
        trmats(ii).Atest10 = out.Atest10;
        trmats(ii)._hptr = _hptr;
     end
    
    fdata.emp10yr   = emp10yr;
    fdata.emp00yr   = emp00yr;
    fdata.f10yr     = f10yr;
    fdata.hpf10yr   = hpf10yr;
    fdata.infocnty  = infocnty;
    fdata.labels    = labels;
    ##fdata.A         = out.A; 
    ##fdata.hptr      = hptr;
    
    #savef = sprintf("%s-%s.txt", savenm, datestr(now, 'yyyymmdd'));
    #printf("saving fdata in %s\n", savenm);
    #save(savef, "fdata", "trmats");
end
