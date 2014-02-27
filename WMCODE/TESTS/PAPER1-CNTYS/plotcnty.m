function out = plotcnty(adata, fips, _lw = 2.5, _prtype = "ppt", _fs=10)
	gray = [0.6 0.6 0.6];

	out = struct();
	cntyi = strmatch(fips, adata.labels); 
	out.cntyi = cntyi;
	out.klerr = adata.klerr(cntyi); 
	out.terr = adata.terr(cntyi); 
	mdata = [adata.emp10yr(1:18,cntyi),  adata.f10yr(1:18,cntyi), adata.emp00yr(1:18,cntyi)];
	fdata = [adata.emp10yr(19:36,cntyi), adata.f10yr(19:36,cntyi), adata.emp00yr(19:36,cntyi)];

	switch _prtype 
        case 'paper'
	        
            figure(1, "visible", "off");
            
	        _s1 = subplot(2,1,1);
	        plot(0:5:85, mdata(:,1), 'k'); # emp
	        hold on;
	        plot(0:5:85, mdata(:,2), 'k-^'); # forecast
	        h = plot(0:5:85, mdata(:,3), 'linewidth', _lw, 'k-o'); # jumpoff
	        set (h, 'color', gray);
	        hold off;
	        _ym = ylim();
	        xlabel ('Age');
	        ylabel ('Population');
	        legend({'Empirical',  'Forecasted', 'Jumpoff'});	
	        title  (sprintf("Male age distribution for FIPS %s", fips));

	        _s2 = subplot(2,1,2);
	        plot(0:5:85, fdata(:,1), 'k');
	        hold on;
	        plot(0:5:85, fdata(:,2), 'k-^');
	        h = plot(0:5:85, fdata(:,3), 'linewidth', _lw, 'k-o');
	        set (h, 'color', gray);
	        hold off;
	        _yf = ylim();
	        xlabel ('Age');
	        ylabel ('Population');
	        legend({'Empirical',  'Forecasted', 'Jumpoff'});
	        title  (sprintf("Female age distribution for FIPS %s", fips));

	        _y = max([_ym; _yf]);
	        ylim(_s1, _y);
	        ylim(_s2, _y);

	        print(sprintf('F%s.png', fips), '-dpng');
            
        case 'ppt'

            _lw = 5.0;
            _fs = 16;
            

            figure(1, "visible", "off");
            
	        _s1 = subplot(2,1,1);
	        plot(0:5:85, mdata(:,1), '-', 'linewidth', _lw);
            _a1 = gca;
            set(_a1, "fontsize", _fs);
	        hold on;
	        plot(0:5:85, mdata(:,2), 'r--', 'linewidth', _lw);
	        h = plot(0:5:85, mdata(:,3), '-', 'linewidth', _lw);
	        set (h, 'color', gray);
	        hold off;
	        _ym = ylim()(2);
		    lh = legend({'Empirical',  'Forecasted', 'Jumpoff'});	
            set(lh, "FontSize", _fs);
	        title  (sprintf("Male FIPS %s", fips));

	        _s2 = subplot(2,1,2);
	        plot(0:5:85, fdata(:,1), '-', 'linewidth', _lw);
            _a2 = gca;
            set(_a2, "fontsize", _fs);
	        hold on;
	        plot(0:5:85, fdata(:,2), 'r--', 'linewidth', _lw);
	        h = plot(0:5:85, fdata(:,3), '-', 'linewidth', _lw);
	        set (h, 'color', gray);
	        hold off;
	        _yf = ylim()(2);
	        lh = legend({'Empirical',  'Forecasted', 'Jumpoff'});
            set(lh, "fontsize", _fs);
	        title  (sprintf("Female FIPS %s", fips));

	        _y = max(_ym, _yf);
	        ylim(_s1, [0 _y]);
	        ylim(_s2, [0 _y]);

            set(_a1, 'ytick', [0 _y/2 _y]);
            set(_a2, 'ytick', [0 _y/2 _y]);
            set(_a1, 'xticklabel', []);

	        print(sprintf('F%s.png', fips), '-dpng');

             
    endswitch
endfunction
