function prob = prob(u,p,sites)global N prob=sigmoid(u(p,sites));     rand('state',sum(100*clock))prob=prob>rand([1 2*N+1]);