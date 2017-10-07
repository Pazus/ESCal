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

    parfor m=1:size(x,4)
        
        x_local = x_short(:,:,:,m);
        y_flipped = y_short(:,:,:,m);

        s3 = size(x_local,3);
        
        Rm = zeros(size(x_local));
        x_reshaped = reshape(x_local,N,[]);
        y_reshaped = reshape(permute(y_flipped,[1,3,2]),[],N);
        
        for i=1:s3
            if strcmpi(shape,'right') %useEnergy
                TempConv = x_reshaped(:,(i-1)*N+1:end)*y_reshaped(1:end-(i-1)*N,:);
            else
                TempConv = x_reshaped(:,1:i*N)*y_reshaped(end-i*N+1:end,:);
            end
            Rm(:,:,i) = TempConv;
        end
        
        if strcmpi(shape,'right') %useEnergy
            R(:,:,:,m) = Rm;
        else
            R(:,:,:,m) = Rm;
        end
    end
    
    if K>0
        if strcmpi(shape,'right') %useEnergy
            R = cat(3,R,zeros(N,N,K,size(x,4)));
        else
            R = cat(3,zeros(N,N,K,size(x,4)),R);
        end
    end
    
end