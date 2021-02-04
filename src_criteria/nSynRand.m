function [nSyn] = nSynRand(VAFCurve,VAFRand)

ths=mean(diff(VAFRand));

nSyn=min(find(diff(VAFCurve)<=0.75*ths));

if isempty(nSyn)
    nSyn=length(VAFCurve);
end