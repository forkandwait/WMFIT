[fdata, trmats, cpop, _cnty ] = runuscntys_2sex_all(STARTIN=3, ENDIN=9, ENDF=4,
                                                     savenm="countyfcst", 
                                                     cntyf="../../DATA/USCOUNTIES/allcntysdata.txt");

[analyze_out err_rank] = analyzeus(fdata, subsetf="none");

plotpaper(analyze_out);
