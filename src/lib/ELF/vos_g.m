function eps=vos_g(z)

    zplus1 = z + 1.0;
    zminus1 = z - 1.0;
    
    dummy1 = log( abs(zplus1) / abs(zminus1));
    dummy2 = angle(zplus1) - angle(zminus1);
    
    outreal = real(z) + 0.5*(1.0 - (real(z)*real(z) - imag(z)*imag(z)))*dummy1 + real(z)*imag(z)*dummy2;
    outimag = imag(z) + 0.5*(1.0 - (real(z)*real(z) - imag(z)*imag(z)))*dummy2 - real(z)*imag(z)*dummy1;

    eps=outreal + 1j*outimag;
end
