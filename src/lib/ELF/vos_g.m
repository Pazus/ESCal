function eps=vos_g(z)

    zplus1 = z + ones(size(z));
    zminus1 = z - ones(size(z));
    
    dummy1 = log( abs(zplus1) ./ abs(zminus1) );
    dummy2 = bsxfun(@minus,angle(zplus1),angle(zminus1));
    
    real_sq = bsxfun(@power,real(z),2);
    imag_sq = bsxfun(@power,imag(z),2);
    
    reim1 = ones(size(z)) - bsxfun(@minus,real_sq,imag_sq);
        
    outreal_1 = real(z) + 0.5*bsxfun(@times,reim1,dummy1);
    outreal = bsxfun(@plus,outreal_1,bsxfun(@times,bsxfun(@times,real(z),imag(z)),dummy2));
    
    outimag_1 = imag(z) + 0.5*bsxfun(@times,reim1,dummy2);
    outimag = bsxfun(@minus,outimag_1,bsxfun(@times,bsxfun(@times,real(z),imag(z)),dummy1));
        
%     outreal = real(z) + 0.5* (1.0 - (real(z)*real(z) - imag(z)*imag(z))) *dummy1 + real(z)*imag(z)*dummy2;
%     outimag = imag(z) + 0.5*(1.0 - (real(z)*real(z) - imag(z)*imag(z)))*dummy2 - real(z)*imag(z)*dummy1;

    eps=outreal+1j*outimag;
end
