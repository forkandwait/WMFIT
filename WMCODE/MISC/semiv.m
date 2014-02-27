function out = semiv (dmat, z, lagw, maxr, varargin)
    ## calculate simple semivariance along with spherical (and other) fits
    ## 
    ## out = semiv (xy, z, lagw, maxr, varargin)

    ## save input
    out.z      = z;
    out.lagw   = lagw;
    out.maxr   = maxr;     
    out.nargin = varargin;

    if find(isnan(out.z))
       error("semiv: nans in z values");
    endif

    ## distance matrix checks and saves
    dmax = max(dmat(:));
    if maxr > dmax
       maxr = ceil(dmax);
       warning("semiv: max radius unnecessarily large, adjusting to %i.", maxr);
    endif 
    out.dmat = dmat;

    ## calculate lagb(oundaries) for lags lagw(idth) maxr(adius)
    out.lagb = (0:lagw:(maxr+lagw))';

    ## set up vector to hold variance for each lag
    out.lagnv   = zeros(length(out.lagb)-1,1);
    out.lagv    = zeros(length(out.lagb)-1,1);
    out.saslagt = zeros(length(out.lagb)-1,1);

    ## for each lag
    for lagi = 1:(length(out.lagb)-1)

        out.saslagt(lagi) = 0;

        ## for each point
        for pnti = 1:rows(dmat) 

            ## for each other point in lag window
            ##     for each index in the distance matrix inside lag
            ## COULD BE VECTORIZED
            lagmatches = (find(dmat(:,pnti) > out.lagb(lagi) & dmat(:,pnti) <= out.lagb(lagi+1)));
            for lagpnti = lagmatches'

                ## collect variance components
                out.lagnv(lagi)  += 1; 
                out.lagv(lagi)   += (z(lagpnti) - z(pnti))^2;    

                ## "the total [not average] distance from the origin Distance=0 of all pairs in a given lag class" ... divide later
                out.saslagt(lagi) += dmat(pnti, lagpnti);

            endfor
        endfor
    endfor
    
    ## delete outermost band if empty
    if out.lagnv(end) == 0 
       out.lagnv(end) = [];
       out.lagv(end)  = []; 
    endif

    ## finish calculation of variance for each lag
    out.lagvt  = out.lagv ./ (2*out.lagnv); 

    ## slightly modified boundaries for plotting
    out.lagb = out.lagb(2:end);

    ## average the total distance from above
    out.saslagb = out.saslagt ./ out.lagnv;
    out.saslagb = out.saslagb;

    ## fit spherical, linear, exponential models to variance
    out.sphf = mksph(out.saslagb, out.lagvt);
    out.optseed = [quantile(out.lagvt, .8), median(out.saslagb)]; # [7 60]
    [out.sphparm, out.fvec, out.info, out.output, out.grad, out.hess] = fminunc(out.sphf{2}, out.optseed);
    if out.info ~= 1
       warning ("semiv: fminunc did not converge. Check starting point, but maybe just best possible. Info code: %i", out.info);
    endif
    out.sphfit = out.sphf{1}(out.sphparm); 
    if 1
        plot(out.saslagb, out.sphfit);
        hold on;
        plot(out.saslagb, out.lagvt, 'ko', 'linewidth', 2.5);
        hold off;
    endif

endfunction

## closure to make the spherical estimate thing from given empirical thing c=x(1) a=x(2)
function outf = mksph(lagb, empv)
    outf{1} = @(x) (ifelse(lagb<=x(2), (x(1)*(1.5*(lagb/x(2)) - 0.5*(lagb/x(2)).^3)), repmat(x(1),rows(lagb),1)));  
    outf{2} = @(x) (sumsq(outf{1}(x) - empv));
endfunction


%!test
%!    pkg load linear-algebra
%!    xy = cartprod(0:2:4, 0:2:4);
%!    dxy = squareform(pdist(xy));
%!    foo = semiv(dxy, rand(9,1), 2, 5);

#{
# coal seam data
coal = [
 0.7  59.6  34.1 
 4.8  52.8  34.3 
 6.4  33.7  36.4  
13.3   0.6  44.7 
17.8   6.9  43.9 
23.0  93.9  43.6 
24.8  26.3  39.7 
27.7  83.3  41.8 
29.5  89.4  43.0 
32.7  40.2  37.5 
37.0  70.3  39.2 
39.4  82.5  41.4 
46.4  84.1  41.5 
51.0  88.8  42.0 
55.5  92.9  42.2 
62.1  26.6  40.1 
70.5  83.7  40.9 
78.1  45.5  38.7 
80.5  55.9  38.7 
84.5  11.0  41.5 
86.7  70.4  39.6 
88.4  12.1  41.3 
88.9   6.2  41.5 
91.5  55.4  39.0 
55.8  50.5  38.1 
 2.1  82.7  42.2 
 5.9  67.1  37.0 
 7.0  46.7  34.6 
13.3  68.2  37.8 
20.1  66.3  37.7 
24.3  73.0  39.3 
26.4  58.0  36.9 
27.9  90.8  43.3 
30.1   6.1  43.6 
34.8   8.1  43.3 
38.2  77.9  40.7 
43.0   4.7  43.3 
46.7  10.6  42.6 
52.8  68.9  39.3 
56.0   1.6  42.7 
63.0  12.7  41.8 
70.9  11.0  41.7 
78.2   9.1  41.7 
81.1  51.0  38.6 
85.2  67.3  39.4 
87.2  55.7  38.8 
88.4  99.6  41.2 
90.6   7.0  41.5 
92.9  46.8  39.1 
96.2  84.3  40.3 
 4.7  75.1  39.5   
 6.0  35.7  35.9   
 8.2  40.1  35.4   
13.4  31.3  37.8   
22.7  87.6  42.8   
24.8  15.1  42.3   
26.9  65.0  37.8   
29.1  47.9  36.7   
30.8  12.1  42.8   
35.3  32.0  38.8   
38.9  23.3  40.5   
43.7   7.6  43.1   
49.9  22.1  40.7   
52.9  32.7  39.2   
60.6  75.2  40.1   
69.0  75.6  40.1   
71.5  29.5  39.8   
78.4  20.0  40.8   
83.8   7.9  41.6   
85.5  73.0  39.8   
88.1   0.0  41.6   
88.8  82.9  40.5   
90.7  49.6  38.9   
93.4  70.9  39.7   
98.2  58.2  39.5   
];
dmat = squareform(pdist(coal(:,1:2)));
semivout = semiv(dmat, coal(:,3), 7, 70);
figure; plot(semivout.saslagb, semivout.lagvt, 'o');


#}
