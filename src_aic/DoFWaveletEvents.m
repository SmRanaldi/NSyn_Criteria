function [out]=DoFWaveletEvents(signalIn,idxEvents)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              %
%   VERSION 2.0 February 2021  %
%                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for cycle=1:length(idxEvents)-1

    signal=signalIn(idxEvents(cycle)+1:idxEvents(cycle+1));
    signal=signal-mean(signal);
    n=floor(log2(length(signal)));
    
    [CA{1},CD{1}]=dwt(signal,'db5');
    
    for i=2:n
        [CA{i},CD{i}]=dwt(CA{i-1},'db5');
    end
    
    rCA=zeros(length(signal),n);
    rCD=zeros(length(signal),n);
    pCA = [];
    pCD = [];

    for i=1:n
        
        rCA(:,i)=resample(CA{i},length(signal),length(CA{i}));
        rCA(:,i)=rCA(:,i)-mean(rCA(:,i));
        rCD(:,i)=resample(CD{i},length(signal),length(CD{i}));
        rCD(:,i)=rCD(:,i)-mean(rCD(:,i));
        
        pCA(:,i)=periodogram(rCA(:,i),[],length(rCA(:,i)));
        pCD(:,i)=periodogram(rCD(:,i),[],length(rCD(:,i)));
        
        pW(1,i)=sum(pCA(:,i));
        pW(2,i)=sum(pCD(:,i));
        
    end
    
    pS=periodogram(signal,[],length(signal));
    
    for i=1:n
        L=round(length(pS)/(2^(i-1)));
        cpCD(i)=corr(pCD(1:L,i)/norm(pCD(1:L,i)),pS(1:L)/norm(pS(1:L)));
    end
    
    [idxCD]=min(find(diff(cpCD)>-0.05*max(abs(diff(cpCD)))));
    if ~isempty(idxCD)
        idx=idxCD(1);        
        if idx==length(CA)            
            idx=idx-1;            
        end        
    else
        idx=length(CA)-1;
    end
    
    out(cycle)=length(CA{idx+1})*(1+min(1,pW(2,idx+1)/pW(1,idx+1)));
    
end