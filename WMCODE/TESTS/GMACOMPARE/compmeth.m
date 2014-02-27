function out = compmeth (fips, P)         

    figure(3); out.outt = singlecnty(fips, mkconwm("tightm"), P, 3, 2);
    figure(2); out.outw = singlecnty(fips, mkconwm("wideopen"), P, 3, 2);
    figure(1); out.outg = singlecnty(fips, mkconwm("generic"), P, 3, 2);
    figure(4); out.outc = singlecnty(fips, mkconwm("college"), P, 3, 2);

    out.mapes = sprintf("g:%f, t:%f, w:%f, c:%f\n", 
                        out.outg.ERR.MedAPE, out.outt.ERR.MedAPE, out.outw.ERR.MedAPE, out.outc.ERR.MedAPE);
    out.kls   = sprintf("g:%f, t:%f, w:%f, c:%f\n", 
                        out.outg.ERR.KLD, out.outt.ERR.KLD, out.outw.ERR.KLD, out.outc.ERR.KLD);
    disp(out.mapes);
    disp(out.kls);

endfunction;


%!demo
%!    P = Pop;
%!    out = compmeth("53001", P);
