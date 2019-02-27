function R = conv3d(x,y,shape, K)
    
    if nargin <3 || isempty(shape)
        shape = 'right';
    end
    
    if nargin <4 || isempty(K)
        K=0;
    end 

    
    N = size(x,1);
        
    if strcmpi(shape,'right') %useEnergy
        x_short = x(:,:,1:end-K,:);
        y_short = flip(y(:,:,1:end-K,:),3);
    else
        x_short = x(:,:,1+K:end,:);
        y_short = flip(y(:,:,1+K:end,:),3);
    end
    
    R = zeros(size(x_short));
    
    s3 = size(x_short,3);    
    M = size(x,4);
    
    x_reshaped = reshape(x_short,N,[],M);
    y_reshaped = reshape(permute(y_short,[1,3,2,4]),[],N,M);
    
    if M>1
        parfor m=1:M
            x_local2 = x_reshaped(:,:,m);
            y_local2 = y_reshaped(:,:,m);            
            R(:,:,:,m) = doCalc(x_local2, y_local2, N, s3, shape);
        end
    else
        R = doCalc(x_reshaped, y_reshaped, N, s3, shape);
    end
    
    if K>0
        if strcmpi(shape,'right') %useEnergy
            R = cat(3,R,zeros(N,N,K,size(x,4)));
        else
            R = cat(3,zeros(N,N,K,size(x,4)),R);
        end
    end
    
end

function Rm = doCalc(x_local2, y_local2, N, s3, shape)
    Rm = zeros(N,N,s3);
    if strcmpi(shape,'right') %useEnergy
        for i=1:s3
            Rm(:,:,i) = x_local2(:,(i-1)*N+1:end)*y_local2(1:end-(i-1)*N,:);
        end
    else
        for i=1:s3
            Rm(:,:,i) = x_local2(:,1:i*N)*y_local2(end-i*N+1:end,:);
        end
    end
end