function H = H(x,t)
global Diff td uc tR 


tt=0.000000000000001;
H = A(x,tt)-A(x,t);


