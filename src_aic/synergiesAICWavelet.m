function [W,H,meanH,VAF,AIC,R2,VAFGlobal,R2Global,DoFTot,VAFMuscles,DoF,L,nPower,E]=synergiesAICWavelet(M,init,criterion,events)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              %
%   VERSION 2.0 February 2021  %
%                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Matrix organization

if size(M,1)>size(M,2)
    M=M';
end
nMuscles=size(M,1);
DoFTot=[];
opt = statset('MaxIter',1000,'TolFun',1e-4,'TolX',1e-4);

%% Factorization and approximation

for i=1:nMuscles
    switch init
        case 'sparse'
            [WProv{i},HProv{i}] = NN_mat_fact_sparse(M,i,100000);
        case 'rand'
            [WProv{i},HProv{i}] = nnmf(M,i,'algorithm','mult','options',opt,'replicates',10);       
    end
    rec{i} = WProv{i}*HProv{i};
    
    %% Signal dependent noise power
    
    MNoise=noiseHVar(HProv{i},WProv{i});
   
    %% Degrees of freedom evaluation
    
    for j=1:i
        [DoFH(j,:)]=DoFWaveletEvents(HProv{i}(j,:),events);
    end
    DoFTot=[DoFTot;DoFH(:)];
    
    %% Model selection criteria
    
    VAFGlobal(i) = 1 - sum(sum(((rec{i} - M).^2)))./sum(sum(M.^2)); % Global Variance Accounted For
    likelihood=[];
    for j=1:nMuscles
        VAFMuscles(i,j)=1-sum((rec{i}(j,:)-M(j,:)).^2)./sum(M(j,:).^2); % Muscle by muscle Variance Accounted For
        E(i,j)=sum((rec{i}(j,:)-M(j,:)).^2);
        likelihood(j)=sum(((rec{i}(j,:) - M(j,:)).^2)./(var(M(:))+MNoise(j,:)));
        nPower(i,j)=mean(MNoise(j,:)); % SDN power
    end
    R2Global(i) = 1 - sum(sum(((rec{i} - M).^2)))./sum(sum((M-mean(mean(M))).^2)); % Global R-squared
    DoF(i)=2*i*(size(M,1))+2*(sum(sum(DoFH))); % Degrees of freedom
    L(i)=sum(likelihood); % Likelihood value
end

%% DoF fit and AIC curve

% [a,b]=ls_fit(linspace(1,nMuscles,nMuscles),DoF); % DoF fit
% AIC=L+a+b.*linspace(1,nMuscles,nMuscles); % AIC value
AIC = L + DoF;

%% Model order selection

nSyn=findNSynAIC(AIC,criterion);

W=WProv{nSyn};
H=HProv{nSyn};
VAF=VAFGlobal(nSyn);
R2=R2Global(nSyn);

%% W and H normalization

for i=1:nSyn
    H(i,:)=H(i,:).*norm(W(:,i));
    W(:,i)=W(:,i)./norm(W(:,i));
end

%% Mean H profile

meanH = []; 
% for i=1:nSamplesCycle
%     meanH(:,i)=mean(H(:,i:nSamplesCycle:end),2);
% end