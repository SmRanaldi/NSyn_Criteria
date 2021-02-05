function [filterCoeff,spectrum] = fdeluca(Fh,Fl,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% function filter = fdeluca(Fh,Fl,fs)
%
%	INPUT	Fh:upper cutoff
%			Fl:lower cutoff
%		    fs: sampling frequency
%
%	OUTPUT
%			filter: coefficients of the filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f = 1:fs/2;
P = (Fh.^4*f.^2)./((f.^2+Fl.^2).*(f.^2+Fh.^2).^2);
M=fs/2;

%calculates the antitransform, starting from the modulus of the spectrum

Px = [P P(M) fliplr(P)];
Px = sqrt(Px);
Px1= Px.*exp(-1i*imag(hilbert(log(Px))));

filterCoeff=real(ifft(Px1,fs));
spectrum = P./length(P);



