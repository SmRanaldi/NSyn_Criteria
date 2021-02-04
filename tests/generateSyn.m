function [M,W,HOut,lengthCycleOut] = generateSyn(nMuscles, nSyn, nSamplesActivation, lengthCycle, nResample)

maxIter=10;

NN = ceil(nMuscles/(nSyn/2));
W = 0.05.*rand(nMuscles,nSyn);
k = randperm(nMuscles);
N=randi([1 NN]);
W(k(1:N),1) = W(k(1:N),1)+0.6+0.05.*rand(N,1);

i=1;iter=0;

while i<nSyn
    
    iter=iter+1;
    
    k = randperm(nMuscles);
    N=randi([1 NN]);
    W(:,i+1)=0.05.*rand(nMuscles,1);
    W(k(1:N),i+1) = W(k(1:N),i+1)+0.6+0.05.*rand(N,1);
    for j=1:i
        ccc(j)=dot(W(:,i+1)./norm(W(:,i+1)),W(:,j)./norm(W(:,j)));
    end
    
    if max(ccc)<0.6
        i=i+1;
        iter=0;
    else if iter==maxIter
            i=i+1;
            iter=0;
        end
    end
    
end

tW=max(W,[],2);

idx=find(tW<0.3);
W(idx,nSyn)=W(idx,nSyn)+0.6*rand(length(idx),1);


H=zeros(nSyn,nSamplesActivation*nSyn);

% for i=1:nSyn
% act=hanning(randi([round(nSamplesCycle/6.5),round(nSamplesCycle/1.5)]))';
% H(i,1:length(act))=act;
% H(i,:)=circshift(H(i,:),sign(randn(1,1))*randi([round(nSamplesCycle/3),10*nSamplesCycle]));
% end
l=nSamplesActivation;
x=1;
sh=randi([-round(l/3),round(l/3)]);
lsh=0;

for i=1:nSyn-1
    a=randi([0,round(nSamplesActivation/2)])+round(nSamplesActivation/4);
    H(i,x-lsh:x-lsh+a-1)=hanning(a)';
    x=x+a-lsh;
end

a=randi([round(nSamplesActivation/2),nSamplesActivation]);
H(nSyn,x-lsh:x-lsh+a-1)=hanning(a)';
H=H(:,1:x-lsh+a-1);

for i=1:nSyn
    H(i,:)=circshift(H(i,:),sh);
end

% for i=1:nSyn
%     H(i,:)=H(i,:).*norm(W(:,i));
%     W(:,i)=W(:,i)./norm(W(:,i));
% end

H=H+0.001*max(H(:));
for i=1:nSyn
    HOut(i,:)=resample(H(i,:),lengthCycle,length(H(i,:)));
    H(i,:)=smooth(H(i,:),ceil(3*(lengthCycle/nResample)));
end

lengthCycleOut=size(HOut,2);

M = W*HOut;