function out = or_compare(orf, or_out);

    for cntyi = 1:length(or_out)
        disp(or_out{cntyi}.fcst.fcstname); 

        if 0
            subplot(2,1,1);
            plot(0:5:85, [orf(1:18,cntyi), or_out{cntyi}.xfcst.Fcst(1:18,2), or_out{cntyi}.fcst.pop(1:18,end), or_out{cntyi}.fcst.pop(1:18,end-2)]);
            legend('or', 'webb', 'true', 'jump');
            
            subplot(2,1,2);
            plot(0:5:85, [orf(19:36,cntyi), or_out{cntyi}.xfcst.Fcst(19:36,2), or_out{cntyi}.fcst.pop(19:36,end), or_out{cntyi}.fcst.pop(19:36,end-2)]);
            legend('or', 'webb', 'true', 'jump');

            input('sdkfjs');
        endif
        pop2010(:,cntyi)  = or_out{cntyi}.fcst.pop(:,end);
        webb2010(:,cntyi) = or_out{cntyi}.xfcst.Fcst(:,2); 
        
    endfor

    out.orf_perr = (pop2010-orf) ./ orf;
    out.webb_perr = (webb2010-orf) ./ orf;
    
    figure(2, 'Visible', true);
    subplot(2,1,1);
    hist(abs(out.orf_perr(:)), 50, 1);
    ylim([0 .25]); xlim([0 0.5]);
    subplot(2,1,2);
    hist(abs(out.webb_perr(:)), 50, 1);
    ylim([0 .25]); xlim([0 0.5]);

    disp('Percent error quantiles');
    qs = [0 .05 .5 .95 1];
    disp([qs', quantile([out.orf_perr(:), out.webb_perr(:)], qs)]);
    disp('Absolute percent error quantiles');
    disp([qs', quantile(abs([out.orf_perr(:), out.webb_perr(:)]), qs)]);
    disp('MAPE');
    disp(mean(abs([out.orf_perr(:), out.webb_perr(:)])));
    disp('MALPE');
    disp(mean(([out.orf_perr(:), out.webb_perr(:)])));

endfunction
