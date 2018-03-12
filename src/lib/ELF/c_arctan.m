function res = c_arctan(z)

	x = real(z);
	y = imag(z);
	imres = imag( -1.0 ./ 4.0 * log(bsxfun(@times,ones(size(x)) - x.^2 - y.^2,ones(size(x)) - x.^2 - y.^2) + 4.0*x.^2) + 1.0 / 2.0 * log(bsxfun(@times,ones(size(y))+y,ones(size(y)) + y) + x.^2));
	if bsxfun(@gt,x,0.0)  %#ok<ALIGN>
		reres = real( repmat(pi/4.0,size(x)) - 0.5 * atan((ones(size(x)) - x.^2 - y.^2) ./ (2.0 * x)));
    elseif bsxfun(@lt,x,0.0)
		reres = real( repmat(pi/4.0,size(x)) - 0.5 * atan((ones(size(x)) - x.^2 - y.^2) ./ (2.0 * x)));
	else
	    reres = repmat( real( pi / 2.0),size(x));	
    end
    res = reres + 1j*imres;
end