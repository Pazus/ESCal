function diimfp = ndiimfp(osc,E0,decdigs,varargin)
%%
%{
   Calculates the normalised DIIMFP
   for a given energy.
%}

if nargin<3, decdigs=10; end
if nargin<2
    warning ('Error in input format')
else
    
%     qmin = sqrt(2*E0/h2ev)-sqrt(2*(E0/h2ev-osc.eloss/h2ev));
%     qmax = sqrt(2*E0/h2ev)+sqrt(2*(E0/h2ev-osc.eloss/h2ev));

    ind = bsxfun(@gt,osc.eloss,0);
    eloss = osc.eloss;    
    qmin = log(sqrt(2*E0/h2ev)-sqrt(2*(E0/h2ev-osc.eloss(ind)/h2ev)));
    qmax = log(sqrt(2*E0/h2ev)+sqrt(2*(E0/h2ev-osc.eloss(ind)/h2ev)));
    q = zeros(length(osc.eloss(ind)),2^(decdigs-1)+1);

%     q = zeros(length(osc.eloss),2^(decdigs-1)+1);
    x_in = zeros(size(osc.eloss));
    
    for i = 1:2^(decdigs-1)+1
        q(:,i) = qmin + (i-1)*(qmax-qmin)/2.^(decdigs-1);
    end
    
    if strcmp( osc.model,'Mermin')
        q(q==0) = 0.01;
    end
    
%     osc.qtran = q/a0;
    osc.qtran = exp(q)/a0;
    osc.eloss = eloss(ind);
    
%     ELF = eps_sum_allwq(osc,'bulk');
%     res = ELF./q;
%     res(isnan(res))=0;
    res = eps_sum_allwq(osc,'bulk');
    
    x_in(1) = 0.0;
    for i=2:length(osc.eloss)
        x_in(i) = 1/pi/(E0/h2ev) * trapz(q(i,:),res(i,:)) *(1/h2ev/a0);
    end

    %% Plot
    %{
    figure;
    xlim([0 100])
    hold on
    box on
    plot(osc.eloss,x_in)
   
    Y = ['biimfp = ',num2str(trapz(osc.eloss,x_in))];
    disp(Y);
    %}

    diimfp = x_in./trapz(eloss,x_in); %normalized
    %diimfp = x_in; %for imfp calculation

end










