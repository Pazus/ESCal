function res = conv_my(x1,x2,dE,shape,indE0)
if nargin<3 || isempty(dE); dE = 1; end;

if ~isscalar(dE); error('Ўаг должен быть скал€ром'); end;

N = numel(x1);
if nargin<5 || isempty(indE0)
    ind = 1:N;
    indE0 = max(ind(x1>0));
    if isempty(indE0); indE0 = N; end;
end
    
if nargin<4 || strcmp(shape,'right')
    
    a = conv(x1,x2);
    res = a(end-(N-indE0)-N+1:end-(N-indE0));
elseif strcmp(shape,'left')
    
    a = conv(x1,x2);
    %res = a(1:numel(x1));
    res = a(1+(N-indE0):N+(N-indE0));
elseif any(strcmp(shape,{'same','full','valid'}))
    res = conv(x1,x2,shape);
end;

res = res *dE;

end