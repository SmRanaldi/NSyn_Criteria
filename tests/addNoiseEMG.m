function y=addNoiseEMG(x,SNR,w,noiseSeries)

idx=find(w>0.001);
Px=var(x(idx));

Pn=Px/(10^(SNR/10));
n=sqrt(Pn)*noiseSeries;
y=x+n;