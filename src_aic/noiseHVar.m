function [out]=noiseHVar(H,W)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              %
%   VERSION 2.0 February 2021  %
%                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(H,1)>size(H,2)
    H=H';
end

for i=1:size(H,1)
    
    H(i,:)=H(i,:).*norm(W(:,i));
    W(:,i)=W(:,i)./norm(W(:,i));
    
end

for i=1:size(H,1)
    
    s=W(:,i)*H(i,:);
        
    nH(i)=std(s(:));
    
end

out=ones(size(H));

for i=1:size(H,1)
    
    out(i,:)=out(i,:)*nH(i);
    
end

out=(W*out).^2;