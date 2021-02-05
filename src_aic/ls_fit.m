function [A,B,r,pCorr] = ls_fit(x,y,varargin)

%%%%%%%%%%%%%%%%%%%%%%%
%[A,B,r,pCorr] = ls_fit(x,y,varargin)
%
% y=A+B*x
%%%%%%%%%%%%%%%%%%%%%%%

nVarargs = length(varargin);
if nVarargs==0
    varplot=0;
else
    varplot=1;
end

if size(x,1)~=size(y,1)
    
    x=x';
    
end

[idxOutliers] = isoutlier(y);
x(idxOutliers) = [];
y(idxOutliers) = [];

A=(sum(x.^2)*sum(y) - sum(x)*sum(x.*y))/(length(x)*sum(x.^2)-(sum(x))^2);
B=(length(x)*sum(x.*y) - sum(x)*sum(y))/(length(x)*sum(x.^2)-(sum(x))^2);

y_est=A+B*x;

[aa,p]=corrcoef(x,y);
r = aa(1,2);
pCorr = p(1,2);
%disp(['R squared: ',num2str(r)]);

if varplot
    scatter(x,y);
    hold on
    plot(x,y_est,'r');
end