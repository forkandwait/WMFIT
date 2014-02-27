function plotpaper(analyze_out);
		 
    close all;
    _g = graphics_toolkit('gnuplot');
    
	fips = {"53009" "53011" "53033" "53063" "53075"  "53077"};
    disp(fips);
	for f = fips
        disp(f{:});
		plotcnty(analyze_out, f{:}, 2.5, "paper");
	endfor
    
    graphics_toolkit(_g);
end
