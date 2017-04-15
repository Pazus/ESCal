function [Xp, Xm]=X_mp(DECS,mu,x,w,M,norm_val)
%X_mp - расчет индикатрисы упругого рассеяния в узлах сетки с разложением
%по азимутальным гармоникам
%   DECS - сечения упругого рассеяния
%   mu   - сетка сечения
%   x    - расчетная сетка задачи, для которой нужно расчитать индикатрису
%   w    - веса расчетной сетки задачи
%   M    - число азимутальных гармоник (если не передается, то ограничиваемся
%   norm_val - нормировочная площадь под индикатрисой
%   нулевой гармоникой)

%% Обработка выходных данных
% if ~isrow(DECS); DECS = DECS.'; end;
% if ~isrow(mu); mu = mu.'; end;
% if ~isrow(x); x = x.'; end;
% if ~isrow(w); w = w.'; end;

%% Расчет сеток
Nk=2.3*1.5;
Lk=1.415; % 1.415;

N = numel(x);                       % Число узлов выходной сетки
N_=floor(Nk*N);                     % Число узлов на временной сетке для получения коэффециентов разложения
L_=floor(Lk*N_); l_=0:L_;           % Число полиномов лежандра

[x_,w_]=legzo_n1(N_);               % Расчет узлов и весов полиномов Лежандра для временной сетки

if nargin<6
    norm_val = 2;
end
if nargin<5
    M = 0;                          % Число азимутальных гармоник по умолчанию
end

if max(mu)<=pi && min(mu)>=0
    mu = cos(mu);
elseif max(mu)>pi && max(mu) <= 180 && min(mu)>=0
    mu = cosd(mu);
elseif min(mu)<-1 || max(mu)>1
    error(['Parameter mu containes unexpected values. min(mu) = ' num2str(min(mu)) ', max(mu) = ' num2str(max(mu))]);
end

%% Инициализация памяти
Xp=zeros(N,N,M+1);                  
Xm=zeros(N,N,M+1);

%% Расчет сечения на временной сетке
if ~isequal(mu,x_)
    DECS = interp1(mu,DECS,x_);
end
%% Расчет полиномов лежандра для временной сетки нулевой гармоники
% PP = Legendre_mu(x_,0,L_);
% xk=(DECS.*w_)*PP';
% xk=norm_val*xk/xk(1);

norm_leg=(2*l_+1)/2;

xk = ExpandFunctionByLegandrePolinomials(x_,DECS,L_,w_);
xk_p = (xk.*norm_leg)';
xk_m = (xk.*(-1).^l_.*norm_leg)';  
% xk_m_norm=repmat(xk_m,1,N);
% xk_p_norm=repmat(xk_p,1,N);
for m=0:M
    P1 = Legendre_mu(x,m,L_);   
    Xp(:,:,m+1) = symmetrize(P1'*diag(xk_p(m+1:end))*P1);
	Xm(:,:,m+1) = symmetrize(P1'*diag(xk_m(m+1:end))*P1);
end
Int=(Xp(:,:,1)+Xm(:,:,1))*w';
Xp=Xp*norm_val./repmat(Int,[1,N,M+1]);
Xm=Xm*norm_val./repmat(Int,[1,N,M+1]);
Xp = symmetrize(Xp);
Xm = symmetrize(Xm);
end