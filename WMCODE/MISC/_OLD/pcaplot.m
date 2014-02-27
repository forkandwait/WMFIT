function pcaplot (u, c1=1, c2=2);
    _r = rows(u)/2;
    plot(u(1:_r,c1), u(1:_r,c2), 'o');
    text(u(1:_r,c1), u(1:_r,c2), num2str([1:_r]'));
endfunction


