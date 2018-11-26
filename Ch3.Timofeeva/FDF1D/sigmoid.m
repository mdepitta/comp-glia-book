function sigmoid = sigmoid(u)

global beta uc

sigmoid=(1.0./(1+exp(-beta*(u-uc)))-1.0/(1+exp(-beta*(0-uc))))./(1-1.0/(1+exp(-beta*(0-uc))));
