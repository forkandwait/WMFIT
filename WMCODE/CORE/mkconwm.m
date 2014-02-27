function constr = mkconwm(modeltype="generic", varargin)
    ## mkcon creates a constraint structure for use in wmh

    ## XXX
    ## In fertility etc there are still some places where indexes are
    ## hardcoded or where appending to the end is not controlled
    ## through a standard variable
    ### warning("mkconwm: Still some naked index numbers!");
		 
	## fiddle with options
	[_regp _propp] = parseparams(varargin);
	if ~isempty(_propp)
		_propp = struct(_propp{:});
	endif

    ## birth ratio
	_r = 0.4886;

	## get matrix infrastructure
	mat = mkmat("twosex");

    ## initial struct with generic constraints
    constr = struct(); 

	## SIMPLE UPPER AND LOWER BOUNDS
	## dl -- diag lower and upper
	dlm = [-0.25 # 1,   0
		   -0.50 # 2,   5
		   -0.50 # 3,  10
		   -1.50 # 4,  15 -- kids are mobile
		   -1.50 # 5,  20 -- kids are mobile
		   -3.50 # 6,  25 -- kids are mobile -- esp home from college
		   -1.50 # 7,  30 -- kids are mobile
		   -0.75 # 8,  35
		   -0.50 # 9,  40
		   -0.50 # 10, 45
		   -0.50 # 11, 50
		   -0.50 # 12, 55
		   -0.75 # 13, 60 -- retirees
		   -0.75 # 14, 65 -- retirees
		   -0.75 # 15, 70 -- retirees
		   -0.50 # 16, 75
		   -0.50 # 17, 80
		   -0.65 # 18, 85
		  ];
	dum = [0.25 # 1,   0, 0.25 
		   0.50 # 2,   5, 0.50 
		   0.50 # 3,  10, 0.50 
		   1.50 # 4,  15, 1.50 -- kids are mobile
		   1.50 # 5,  20, 1.50 -- kids are mobile
		   3.50 # 6,  25, 3.50 -- kids are mobile -- esp home from college
		   1.50 # 7,  30, 1.50 -- kids are mobile
		   0.75 # 8,  35, 0.75 
		   0.50 # 9,  40, 0.50 
		   0.50 # 10, 45, 0.50 
		   0.50 # 11, 50, 0.50 
		   0.50 # 12, 55, 0.50 
		   0.75 # 13, 60, 0.75 -- retirees
		   0.75 # 14, 65, 0.75 -- retirees
		   0.75 # 15, 70, 0.75 -- retirees
		   0.50 # 16, 75, 0.50 
		   0.50 # 17, 80, 0.50 
		   0.65 # 18, 85, 0.65 
		  ];
	dlf = [-0.25 # 1,   0
		   -0.50 # 2,   5
		   -0.50 # 3,  10
		   -1.50 # 4,  15 -- kids are mobile
		   -3.50 # 5,  20 -- kids are mobile
		   -3.50 # 6,  25 -- kids are mobile -- esp home from college
		   -1.50 # 7,  30 -- kids are mobile
		   -0.75 # 8,  35
		   -0.50 # 9,  40
		   -0.50 # 10, 45
		   -0.50 # 11, 50
		   -0.50 # 12, 55
		   -0.75 # 13, 60 -- retirees
		   -0.75 # 14, 65 -- retirees
		   -0.75 # 15, 70 -- retirees
		   -0.50 # 16, 75
		   -0.50 # 17, 80
		   -0.65 # 18, 85
		  ];
	duf = [0.25 # 1,   0
		   0.50 # 2,   5
		   0.50 # 3,  10
		   1.50 # 4,  15 -- kids are mobile
		   1.50 # 5,  20 -- kids are mobile
		   3.50 # 6,  25 -- kids are mobile -- esp home from college
		   1.50 # 7,  30 -- kids are mobile
		   0.75 # 8,  35
		   0.50 # 9,  40
		   0.50 # 10, 45
		   0.50 # 11, 50
		   0.50 # 12, 55
		   0.75 # 13, 60 -- retirees
		   0.75 # 14, 65 -- retirees
		   0.75 # 15, 70 -- retirees
		   0.50 # 16, 75
		   0.50 # 17, 80
		   0.65 # 18, 85
		  ];
	## sdl -- subdiag lower and upper (from HMD m70, m07, f70, f07)
	sd = [0.99614   0.99887   0.99709   0.99908        # 1, 0y
		  0.99769   0.99921   0.99854   0.99938        # 2, 5
		  0.99523   0.99775   0.99784   0.99885        # 3, 10
		  0.98988   0.99380   0.99645   0.99781        # 4, 15
		  0.98952   0.99285   0.99605   0.99741        # 5, 20
		  0.98955   0.99283   0.99503   0.99680        # 6, 25
		  0.98676   0.99178   0.99258   0.99561        # 7, 30
		  0.98066   0.98877   0.98875   0.99333        # 8, 35
		  0.96986   0.98278   0.98285   0.98954        # 9, 40
		  0.95364   0.97377   0.97465   0.98425        # 10, 45
		  0.92778   0.96162   0.96300   0.97767        # 11, 50
		  0.89199   0.94611   0.94672   0.96668        # 12, 55
		  0.84307   0.92023   0.92108   0.94846        # 13, 60
		  0.78197   0.88355   0.87881   0.92151        # 14, 65
		  0.69812   0.82576   0.80976   0.87797        # 15, 70
		  0.59611   0.73736   0.71033   0.80798        # 16, 75
		  0.46707   0.61004   0.57199   0.69610];      # 17, 80
	sdlm = sd(:,1);
	sdum = sd(:,2);
	sdlf = sd(:,3);
	sduf = sd(:,4);
	
	## ... setup w/ defaults = 0 ..Inf	
	constr.lb = nan(length(mat.nz), 1);
 	constr.lb(:) = -Inf;
	constr.ub = nan(length(mat.nz), 1);
 	constr.ub(:) = Inf;
    
    ##     ... subd constraints ... 
	constr.lb(mat._sim) = sdlm; 
	constr.lb(mat._sif) = sdlf;
	constr.ub(mat._sim) = sdum;
	constr.ub(mat._sif) = sduf;  
	
    ##     ... diag constraints (age zero special)... 
	constr.lb(mat._dim) = dlm;
	constr.lb(mat._dif) = dlf;
	constr.ub(mat._dim) = dum;
	constr.ub(mat._dif) = duf;


	## EQUALITY CONSTR (use negatives for lower bounds)

	## empty setup to be added to.
	Aeq = zeros(0, length(mat.nz));
	beq = zeros(0, 1);

    ## fertility proportions by age, transpose to fit Aeq
	if 1
		Fpr = [0.05      # 0-5 (out of region moms 1
			   0.03      # 10-15                   2
			   0.20      # 15-20                   3
			   0.27      # 20-25                   4
			   0.25      # 25-30                   5
			   0.15      # 30-35                   6
			   0.04      # 35-40                   7
			   0.01      # 40-45                   8
			   0.00]';   # 45-50                   9
		
		##     ... renormalize so can adjust above safely ...
		if sum(Fpr) ~= 1.0
			warning("mkconwm: Fertilility prop not sum to 1.0, adjusting");
			Fpr = Fpr/sum(Fpr,2); 
		end
		
		##     ... apply sex ratios ...
		fertpr = [(1-_r)*Fpr _r*Fpr];
		
		##     .. assemble and append to equality matrices...
		_Aeq = zeros(0, length(mat.nz));
		_fitp = [mat._dim(1); mat._fim; mat._dif(1); mat._fif];
		_L = length(_fitp);
		for ii = 1:(_L-1)
			_p = fertpr(ii);
			_v = repmat(-_p, _L, 1);
			_v(ii) = 1-_p;
			_Aeq(ii, _fitp) = _v;
		end
		Aeq = [Aeq; _Aeq];
		beq = [beq; zeros(rows(_Aeq),1)]; 
	endif
	
	## Put into return struct
	constr.Aeq = Aeq;
	constr.beq = beq;


	## INEQUALITY CONSTRAINTS

	## Start basic framework
	A_in = zeros(0, length(mat.nz));
	a_ub = zeros(0,1);

	##     ... monotonic survival ...
	if strmatch(modeltype, "generic")
	   for sdi = 2:16
		   _A_in = zeros(1,length(mat.nz));
		   _A_in(mat._sim(sdi)) = -1;
		   _A_in(mat._sim(sdi+1)) = 1;
		   A_in = [A_in; _A_in];
		   a_ub = [a_ub; -0.0001];

		   _A_in = zeros(1,length(mat.nz));
		   _A_in(mat._sif(sdi)) = -1;
		   _A_in(mat._sif(sdi+1)) = 1;
		   A_in = [A_in; _A_in];
		   a_ub = [a_ub; -0.0001];
		endfor
	endif

	##     ... total fert constraint ...
	if 1
		constr.fsumub = 4.00;
		constr.fsumlb = 0.50;
		_A_in = zeros(2,length(mat.nz));
		_A_in(1, [mat._fim', mat._fif']) = 1;
		_A_in(2, [mat._fim', mat._fif']) = -1;
		A_in = [A_in; _A_in];
		a_ub = [a_ub; constr.fsumub];
		a_ub = [a_ub; -constr.fsumlb];
	endif

	## put stuff in the return struct
	constr.A_in = A_in;
	constr.a_ub = a_ub;
    constr.a_lb = [];

	## TWEAK BASED ON TYPE
    switch modeltype
		case "generic"
			 ;

        case "exactish"
            warning("mkconwm: Refactor me! \"exactish\" modeltype rates-to-cells mathematics needs work! Also needs variable indices!");

            ## get mortality and fertility from property value thing.
            
            if isfield(_propp, "Fert")
                ## fertility -- calc spine, calc range +/-, set total
                ## from spine. XXX not working with saved indexes
                ## here.

                _fertp = [0.001; getfield(_propp, "Fert")];
                if ! all(size(_fertp) == [9 1] )
                    error("mkconwm: bad size for fert: %i x %i; should be 8 x 1.", size(getfield(_propp, "Fert")));
                endif 

                
                ## ... fert totals ...
                constr.a_ub(end-1) = sum(_fertp);
                constr.a_ub(end)   = -sum(_fertp);

                ## ... fert proportions -- based on orig non-parameterized
                ## code above -- should be function ...
                Fpr = _fertp ./ sum(_fertp);
		        fertpr = [(1-_r)*Fpr _r*Fpr]; 
		        _Aeq = zeros(0, length(mat.nz));
		        _fitp = [mat._dim(1); mat._fim; mat._dif(1); mat._fif];
		        _L = length(_fitp);
		        for ii = 1:(_L-1)
			        _p = fertpr(ii);
			        _v = repmat(-_p, _L, 1);
			        _v(ii) = 1-_p;
			        _Aeq(ii, _fitp) = _v;
		        end
		        Aeq = _Aeq;
		        beq = zeros(rows(_Aeq),1);

                ## ... save it.  DANGER if ever use Aeq and beq for something besides Fert
	            constr.Aeq = Aeq;
	            constr.beq = beq;
            endif
            
            ## Mortality -- use upper and lower bounds based on input spine + wiggle
            if isfield(_propp, "Mort") 
                _mort = getfield(_propp, "Mort");
                if ! all(size(_mort) == [18 2] )
                    error("mkconwm: bad size for mortality: %i x %i; should be 18 x 2.", size(_mort));
                endif
                _mortwiggle = 0.05;
	            constr.lb(mat._sim) = _mort(:,1) - _mortwiggle; 
	            constr.lb(mat._sif) = _mort(:,2) - _mortwiggle;
	            constr.ub(mat._sim) = _mort(:,1) + _mortwiggle; 
	            constr.ub(mat._sif) = _mort(:,2) + _mortwiggle; 
            endif

		case "college"
			constr.lb(mat._sim(5)) = 0.4;
			constr.lb(mat._sif(5)) = 0.4;
             
        case "wideopen"
            ##     ... subd constraints ... 
	        constr.lb(mat._sim) = -2.0; 
	        constr.lb(mat._sif) = -2.0;
	        constr.ub(mat._sim) = 3.0;
	        constr.ub(mat._sif) = 3.0;

        case "tightm"
            ##     ... d constraints ... 
             _tm = 0.25;
	        constr.lb(mat._dim(1:end-1)) = -_tm; 
	        constr.lb(mat._dif(1:end-1)) = -_tm;
	        constr.ub(mat._dim(1:end-1)) = _tm;
	        constr.ub(mat._dif(1:end-1)) = _tm;

		otherwise
			error("%s: no such modeltype: %s", mfilename, modeltype);
    end


    ## STARTING GUESS
	x0(mat._sim) = mean([constr.lb(mat._sim), constr.ub(mat._sim)], 2);
	x0(mat._sif) = mean([constr.lb(mat._sif), constr.ub(mat._sif)], 2);
	x0(mat._dim) = mean([constr.lb(mat._dim), constr.ub(mat._dim)], 2); 
	x0(mat._dif) = mean([constr.lb(mat._dif), constr.ub(mat._dif)], 2);
	x0([mat._dim(1); mat._fim; mat._dif(1); mat._fif]) = fertpr; # overwrites corner
	x0 = x0';
	constr.x0 = x0;

end;

