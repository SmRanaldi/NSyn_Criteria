function [nSyn]=nSyn5Perc(VAFIn)

dVAF=diff(VAFIn);

nSyn=min(find(dVAF<=0.05));

if isempty(nSyn)
    nSyn=length(VAFIn);
end