function [W,H,M,EMG,events,lCycle,varOut]=genEMGSyn(nMuscles,lengthActivation,nSyn,nCycles,SNRLimits)

W=0.1*rand(nMuscles,nSyn);

act=randperm(nMuscles);

for i=1:nSyn
    W(act(i),i)=W(act(i),i)+0.65+0.2*rand(1);
end

if nSyn<nMuscles
    
    act=act(nSyn+1:end);
    
    for i=1:length(act)
        
        for j=1:randi(2,1,1)
            
            aa=randi(nSyn,1,1);
            
            W(act(i),aa)=0.65+0.2*rand(1);
            
        end
        
    end
    
end

% for i=1:nSyn
%
%     W(randi(nMuscles,1),i)=0.2+0.4*rand(1);
%
% end

L=lengthActivation;
H=[];

% H(1,:)=generateHLine(L,nSyn);

i=1;
it=0;
maxIter=1000;

while i<=nSyn
    
    L=randi(round([lengthActivation/2;(3/2)*lengthActivation]));
    
    HProv=zeros(nSyn,L);
    HProv(i,1:L-round(L/5))=hanning(L-round(L/5));
%     HProv(i,:)=circshift(HProv(i,:),randi([-round(L/6) round(L/6)]));
    
    H=[H,HProv];
    
    %     it=it+1;
    
    %     H(i,:)=generateHLine(L,nSyn);
    
    %     for j=1:i-1
    %         dotPr(j)=dot(H(i,:)/norm(H(i,:)),H(j,:)/norm(H(j,:)));
    %     end
    
    %     if max(dotPr)<0.4 || it>=maxIter
    i=i+1;
    %         it=0;
    %     end
    
    %     dotPr=[];
    
end

HProv=H;
H=[];

% for i=1:nSyn
%     HProv(i,:)=circshift(HProv(i,:),((-1)^(randi([1,2])))*randi([0,round(lengthActivation/8)]));
% end
events=1;
for i=1:nCycles
    
    H=[H,HProv];
    events=[events;size(H,2)];
    
end

for i=1:nSyn %% Every synergy explains the same amount of variance
    
    temp=W(:,i)*H(i,:);
    H(i,:)=H(i,:)./(nSyn*std(temp(:)));
    %
    %     W(:,i)=W(:,i)/norm(W(:,i));
    %     H(i,:)=H(i,:)/sum(H(i,:));
    temp=W(:,i)*H(i,:);
    varOut(i)=var(temp(:));
    
end

M=W*H;

M=M+(( (0.15*rand(size(M))+0.1) .*M).^2).*randn(size(M));

M(M<0.001)=0.001;

SNR=rand(1,nMuscles);
SNR=(SNR.*(max(SNRLimits)-min(SNRLimits)))+min(SNRLimits(1));

for i=1:nMuscles
    EMG(i,:)=simulateEMG(450,20,1000,M(i,:),randn(size(M(i,:))));
    EMG(i,:)=addNoiseEMG(EMG(i,:),SNR(i),M(i,:),randn(size(M(i,:))));
end

lCycle=round(mean(diff(find(events))));