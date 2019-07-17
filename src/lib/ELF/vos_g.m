function [out, outim]=vos_g(z, img_z)

    if nargin < 2
        reZ = real(z);
        imgZ = imag(z);
    else
        reZ = z;
        imgZ = img_z;
    end

    zplus1 = z + 1;
    zminus1 = z - 1;
    
    if any(any(imgZ ~= 0))
        imgZ2 = imgZ.*imgZ;
        % dummy1 = log( abs(zplus1) ./ abs(zminus1) );
        dummy1 = log( sqrt((zplus1.*zplus1+imgZ2) ./ (zminus1.*zminus1+imgZ2)) );
    %     dummy2 = angle(zplus1) - angle(zminus1);
        dummy2 = atan2(imgZ, zplus1) - atan2(imgZ, zminus1);

        reim1 = 1 - (reZ.^2 - imgZ2);

        outreal_1 = reZ + 0.5*reim1.*dummy1;
        outreal = outreal_1 + reZ.*imgZ.*dummy2;

        outimag_1 = imgZ + 0.5*reim1.*dummy2;
        outimag = outimag_1 - reZ.*imgZ.*dummy1;
    else
%         imgZ2 = imgZ.*imgZ;
        dummy1 = log( abs(zplus1) ./ abs(zminus1) );
    %     dummy2 = angle(zplus1) - angle(zminus1);
        dummy2 = atan2(0, zplus1) - atan2(0, zminus1);

        reim1 = 1 - reZ.^2;

        outreal_1 = reZ + 0.5*reim1.*dummy1;
        outreal = outreal_1;%; + reZ.*imgZ.*dummy2;

        outimag_1 = 0.5*reim1.*dummy2;% + imgZ
        outimag = outimag_1;% - reZ.*imgZ.*dummy1;        
    end
        
%     outreal = real(z) + 0.5* (1.0 - (real(z)*real(z) - imag(z)*imag(z))) *dummy1 + real(z)*imag(z)*dummy2;
%     outimag = imag(z) + 0.5* (1.0 - (real(z)*real(z) - imag(z)*imag(z))) *dummy2 - real(z)*imag(z)*dummy1;

    if nargout == 1
        out = complex(outreal,outimag);
    else
        out = outreal;
        outim = outimag;
    end
end
