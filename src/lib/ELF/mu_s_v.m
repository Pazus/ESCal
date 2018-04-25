function [res]=mu_s_v(osc,a,b,decdigs,w_current,depth,alpha,E0,term,varargin)

if nargin<9, decdigs=10; end
if nargin<8
   warning ('Error in input format')
else
   if a==b
       osc.qtran = a;
   else
       osc.qtran = a:(b-a)/2^(decdigs-1):b;
   end
   
   osc.eloss = w_current;
   theta = 0:pi/4:pi/2;
   phi = 0:2*pi/5:2*pi;
   
   switch(term)
       case 1
           romall = eps_sum(osc)./(osc.qtran*a0);
       case 2
           qz = bsxfun(@times,osc.qtran*a0,cos(theta)');
           Q = bsxfun(@times,osc.qtran*a0,sin(theta)');
           v_per = cosd(alpha).*sqrt(2*E0/h2ev);
           r = bsxfun(@rdivide,depth/a0,cosd(alpha));
           
           qzrcos = bsxfun(@times,qz,bsxfun(@times,r,cosd(alpha))); % qz*r*cos(alpha)
           Qcosalpha = bsxfun(@times,Q,cosd(alpha)); % Q*cos(alpha)
           excos = bsxfun(@times,cos(qzrcos),exp(bsxfun(@times,(-1)*abs(r),Qcosalpha)));
           
           weight = bsxfun(@times,bsxfun(@times,osc.qtran*a0,(sin(theta).^2)'),excos);
           qsintheta = bsxfun(@times,(osc.qtran*a0).*sqrt(2*E0/h2ev),sin(theta)');
           phialpha = bsxfun(@times,cos(phi),sind(alpha));
           
           exdim = repmat(qsintheta, 1, 1, length(phialpha)); %add extra dimension over phi
           B = bsxfun(@times,exdim,reshape(phialpha,1,1,[]));
           
           Qv_per = repmat(bsxfun(@times,Q.^2,v_per.^2),1, 1, length(phialpha)); %add extra dimension over phi
           w_wave = bsxfun(@minus,osc.eloss/h2ev,B);
           
           osc.qtran = Q(1,:);
           Im(1,:) = eps_sum(osc);
           osc.qtran = Q(2,:);
           Im(2,:) = eps_sum(osc);
           osc.qtran = a:(b-a)/2^(decdigs-1):b;
           
           del = bsxfun(@plus,w_wave.^2,Qv_per);
           romall = bsxfun(@times,Im,bsxfun(@rdivide,repmat(weight,1, 1, length(phialpha)),del));
       case 3
           Im = zeros(length(theta),length(osc.qtran));
           qz = bsxfun(@times,osc.qtran.*a0,cos(theta)');
           Q = bsxfun(@times,osc.qtran.*a0,sin(theta)');
           v_per = cosd(alpha).*sqrt(2*E0/h2ev);
           r = bsxfun(@rdivide,depth/a0,cosd(alpha));
           
           qzrcos = bsxfun(@times,qz,bsxfun(@times,r,cosd(alpha))); % qz*r*cos(alpha)
           Qcosalpha = bsxfun(@times,Q,cosd(alpha)); % Q*cos(alpha)           
           qsintheta = bsxfun(@times,osc.qtran.*a0,(sin(theta).^2)');
           
           top = bsxfun(@times,bsxfun(@times,qsintheta,cos(qzrcos)),exp(bsxfun(@times,(-1)*abs(r),Qcosalpha)));
           
           qvsintheta = bsxfun(@times,(osc.qtran*a0).*sqrt(2*E0/h2ev),sin(theta)');
           exdim = repmat(qvsintheta, 1, 1, length(phi)); %add extra dimension over phi
           B = bsxfun(@times,exdim,reshape(cos(phi),1,1,[]));
           
           w_wave = bsxfun(@minus,osc.eloss/h2ev,B.*sind(alpha));
           Qv_per = bsxfun(@times,Q.^2,v_per.^2); %add extra dimension over phi
           bottom = bsxfun(@plus,w_wave.^2,Qv_per);

           for i = 1:length(theta)
               osc.qtran = Q(i,:)/a0;
               Im(i,:) = eps_sum1(osc);
           end
           osc.qtran = a:(b-a)/2^(decdigs-1):b;

           romall = bsxfun(@times,Im,bsxfun(@rdivide,top,bottom));
       case 4
           Im = zeros(length(theta),length(osc.qtran));
           
           Q = bsxfun(@times,osc.qtran.*a0,sin(theta)');
           v_per = cosd(alpha).*sqrt(2*E0/h2ev);
           r = bsxfun(@rdivide,depth/a0,cosd(alpha));
           
           Qcosalpha = bsxfun(@times,Q,cosd(alpha)); % Q*cos(alpha)           
           qsintheta = bsxfun(@times,osc.qtran.*a0,(sin(theta).^2)');
           
           top = bsxfun(@times,qsintheta,exp(bsxfun(@times,(-1)*abs(r),Qcosalpha)));
           
           qvsintheta = bsxfun(@times,(osc.qtran*a0).*sqrt(2*E0/h2ev),sin(theta)'); %q*v*sin(theta)
           exdim = repmat(qvsintheta, 1, 1, length(phi)); %add extra dimension over phi
           B = bsxfun(@times,exdim,reshape(cos(phi),1,1,[]));
           
           w_wave = bsxfun(@minus,osc.eloss/h2ev,B.*sind(alpha));
           Qv_per = bsxfun(@times,Q.^2,v_per.^2); %add extra dimension over phi
           bottom = bsxfun(@plus,w_wave.^2,Qv_per);

           for i = 1:length(theta)
               osc.qtran = Q(i,:)/a0;
               Im(i,:) = eps_sum1(osc);
           end
           osc.qtran = a:(b-a)/2^(decdigs-1):b;

           temp = bsxfun(@times,Im,bsxfun(@rdivide,top,bottom));
           
           newexp = exp(bsxfun(@times,(-1)*abs(r),Qcosalpha));
           add = bsxfun(@minus,2.*cos(w_wave.*r./sqrt(2*E0/h2ev)),newexp);
           
           romall = bsxfun(@times,temp,add);
   end
   romall(isnan(romall))=0;   

   switch term
       case 1
           res = trapz(osc.qtran*a0,romall);
       otherwise
           res = trapz(osc.qtran*a0,trapz(theta,trapz(phi,romall,3),1));
   end

end