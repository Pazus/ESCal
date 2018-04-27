function [diimfp,dsep] = ndiimfp_Li(osc,E0,depth,alpha,decdigs,varargin)
%%
%{
   Calculates the normalised NDIIMFP (in eV^-1*A^-1)
   for a given energy, angle and depth
   from solid to vacuum
   according to the Li algorithm Eq.(9)
   Y.C. Li et al. / Surface Science 589 (2005) 67-76.
%}

if nargin<5, decdigs=10; end
if nargin<4
    warning ('Error in input format')
else
    
    qmin = sqrt(2*E0/h2ev)-sqrt(2*(E0/h2ev-osc.eloss/h2ev));
    qmax = sqrt(2*E0/h2ev)+sqrt(2*(E0/h2ev-osc.eloss/h2ev));
    
    q = zeros(length(osc.eloss),2^(decdigs-1)+1);
    
    for i = 1:2^(decdigs-1)+1
        q(:,i) = qmin + (i-1)*(qmax-qmin)/2.^(decdigs-1);
    end
    
    theta = 0:pi/6:pi/2;
    phi = 0:2*pi/10:2*pi;
    
    Im = zeros(length(osc.eloss),2^(decdigs-1)+1,length(theta));
    
    qz = bsxfun(@times,repmat(q, 1, 1, length(theta)),reshape(cos(theta),1,1,[]));
    
    Q = bsxfun(@times,repmat(q, 1, 1, length(theta)),reshape(sin(theta),1,1,[]));
    v_per = cosd(alpha).*sqrt(2*E0/h2ev);
    r = depth/a0/cosd(alpha);
    
    qzrcos = qz.*(r*cosd(alpha)); % qz*r*cos(alpha)
    Qcosalpha = Q.*cosd(alpha); % Q*cos(alpha)
    qsintheta = bsxfun(@times,repmat(q, 1, 1, length(theta)),reshape(sin(theta).^2,1,1,[]));
    
    qvsintheta = bsxfun(@times,repmat(q.*sqrt(2*E0/h2ev), 1, 1, length(theta)),reshape(sin(theta),1,1,[]));
    exdim = repmat(qvsintheta, 1, 1, 1, length(phi)); %add extra dimension over phi
    B = bsxfun(@times,exdim,reshape(cos(phi),1,1,1,[]));
    
    w_wave = bsxfun(@minus,repmat((osc.eloss/h2ev)',1,2^(decdigs-1)+1,length(theta),length(phi)),B.*sind(alpha));
    Qv_per = bsxfun(@times,Q.^2,v_per.^2); %add extra dimension over phi
    bottom = bsxfun(@plus,w_wave.^2,Qv_per);
    
    %% Bulk
    x_in_clear_b = zeros(size(osc.eloss));
    x_in_b = zeros(size(osc.eloss));
    
    if depth<=0
        % clear bulk
        osc.qtran = q/a0;
        ELF=eps_sum_test(osc,'bulk');
        res = ELF./q;
        res(isnan(res))=0;
        for i=1:length(osc.eloss)
            x_in_clear_b(i) = 1/pi/(E0/h2ev)*stepfunction(-depth) * trapz(q(i,:),res(i,:))/h2ev/a0;
        end
        
        % reduced bulk
        top_in = bsxfun(@times,bsxfun(@times,qsintheta,cos(qzrcos)),exp(bsxfun(@times,(-1)*abs(r),Qcosalpha)));
        
        for i = 1:length(theta)
            osc.qtran = Q(:,:,i)/a0;
            Im(:,:,i) = eps_sum_test(osc,'bulk');
        end
        
        romall_in = bsxfun(@times,Im,bsxfun(@rdivide,top_in,bottom));
        romall_in(isnan(romall_in))=0;
        res_in = trapz(theta,trapz(phi,romall_in,4),3);
        
        for i=1:length(osc.eloss)
            x_in_b(i) = (-2)*cosd(alpha)/(pi^3)*stepfunction(-depth) * trapz(q(i,:),res_in(i,:)) /h2ev/a0;
        end
    end
    
    %% Surface
    x_in = zeros(size(osc.eloss));
    res_in = zeros(length(osc.eloss),2^(decdigs-1)+1);
    res_out = zeros(length(osc.eloss),2^(decdigs-1)+1);
    
    for i = 1:length(theta)
        osc.qtran = Q(:,:,i)/a0;
        Im(:,:,i) = eps_sum_test(osc,'surface');
    end
    
    if depth<=0
        top_in = bsxfun(@times,bsxfun(@times,qsintheta,cos(qzrcos)),exp(bsxfun(@times,(-1)*abs(r),Qcosalpha)));
        romall_in = bsxfun(@times,Im,bsxfun(@rdivide,top_in,bottom));
        romall_in(isnan(romall_in))=0;
        res_in = trapz(theta,trapz(phi,romall_in,4),3);
    end
    
    if depth>=0
        top_out = bsxfun(@times,qsintheta,exp(bsxfun(@times,(-1)*abs(r),Qcosalpha)));
        newexp = exp(bsxfun(@times,(-1)*abs(r),Qcosalpha));
        add = bsxfun(@minus,2.*cos(w_wave.*r./sqrt(2*E0/h2ev)),newexp);
        
        romall_out = bsxfun(@times,bsxfun(@times,Im,bsxfun(@rdivide,top_out,bottom)),add);
        romall_out(isnan(romall_out))=0;
        res_out = trapz(theta,trapz(phi,romall_out,4),3);
    end
    
    for i=1:length(osc.eloss)
        x_in(i) = 4*cosd(alpha)/(pi^3)/h2ev/a0 * (trapz(q(i,:),res_in(i,:))* stepfunction(-depth) + trapz(q(i,:),res_out(i,:))* stepfunction(depth));
    end
    
    %% Plot
    figure;
    xlim([0 100])
    hold on
    box on
    plot(osc.eloss,x_in_b + x_in_clear_b)
    plot(osc.eloss,x_in)
    plot(osc.eloss,x_in_clear_b + x_in_b + x_in)
    legend('Reduced bulk','Surface','DIIMFP');
    
    Y = ['siimfp = ',num2str(trapz(osc.eloss,x_in))];
    disp(Y);
    Y = ['biimfp = ',num2str(trapz(osc.eloss,x_in_b + x_in_clear_b))];
    disp(Y);
    
    if depth==0
        diimfp = x_in_clear_b./0.5./trapz(osc.eloss,x_in_clear_b./0.5);
    else
        diimfp = x_in_clear_b./trapz(osc.eloss,x_in_clear_b);
    end
    dsep = x_in./trapz(osc.eloss,x_in);
%     dsep = (x_in_clear_b + x_in_b + x_in)./trapz(osc.eloss,x_in_clear_b + x_in_b + x_in);

end

%% Heaviside function
function x = stepfunction(depth)
    if depth > 0
        x = 1;
    elseif depth < 0
        x = 0;
    elseif depth==0
        x = 0.5;
    end
end

end










