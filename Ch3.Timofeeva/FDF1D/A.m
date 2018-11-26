function A = A(x,t)
global Diff td uc tR

A = sqrt(td*Diff)/(4*Diff)*(exp(-abs(x)/sqrt(td*Diff)).*erfc(-abs(x)./sqrt(4*Diff*t)+sqrt(t/td))+ ...
    exp(abs(x)/sqrt(td*Diff)).*erfc(abs(x)./sqrt(4*Diff*t)+sqrt(t/td))).*heav(t);

