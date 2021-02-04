function [W,emgs,H,M,lengthCycle]=generateEMGSynergies(nMuscles,nSyn,SNRValues,nSamplesActivation,lengthCycle,nCycles,nResample)

[MProv,W,H,lengthCycle]=generateSyn(nMuscles,nSyn,nSamplesActivation,lengthCycle,nResample);
M=[];

SNR=rand(1,nMuscles);
SNR=(SNR.*(max(SNRValues)-min(SNRValues)))+min(SNRValues(1));

for i=1:nCycles
    M=[M,MProv];
end

for i=1:nMuscles
    M(i,:)=M(i,:)./norm(M(i,:));
end

for i=1:nMuscles
    emgs(i,:)=simulateEMG(450,20,1000,M(i,:),randn(size(M(i,:))));
    emgs(i,:)=addNoiseEMG(emgs(i,:),SNR(i),M(i,:),randn(size(M(i,:))));
end
