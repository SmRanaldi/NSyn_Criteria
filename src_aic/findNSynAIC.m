function [N]=findNSynAIC(AIC,criterion)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              %
%   VERSION 2.0 February 2021  %
%                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch criterion
    
    case 'min'
        
        %% Minimum value
        
        [~, N]=min(AIC);
        
    case 'der'
        
        %% First stationary point
        
        dAIC=diff(AIC);
        
        ths(1)=abs(dAIC(1));
        
        for i=2:length(dAIC)
            
            ths(i)=max(abs(dAIC(1:i-1)));
            
        end
                              
        dAIC=dAIC./ths;
        
        iAIC=find(dAIC>-0.005);
        
        N=iAIC(1);
                
    case 'firstpeak'
        
        %% First negative peak
        
        [~, idx]=findpeaks(-AIC);
        
        N=idx(1);
        
    case 'lastpeak'
        
        %% Last negative peak
        
        [~, idx]=findpeaks(-AIC);
        
        N=idx(end);
        
end