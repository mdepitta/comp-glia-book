% One-dimensional case. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear

global Diff td uc tR beta refractory_max N sigma 

d = 2.0*10^(-6);
Diff = 30*10^(-12); % d*d;
td = 0.2;     
tR = 0.01;  
uc = 0.1;
sigma=1;
beta = 100000;  

rand('state',sum(100*clock))

refractory_max = 70;

N = 25;   % number of calcium stores
L = 26; 

s = 10;    % number of mesh points between each calcium store.
step = d/s;   % mesh for x; 
x = linspace(-L*d,L*d,2*L*d/step+1);
st = linspace(-N*d,N*d,2*N+1);

p_max = 100;

for i=1:2*N+1
    sites(i) = (L-N)*d/step + 1 + s*(i-1);
end

a = zeros(p_max,length(sites));

uend = zeros(p_max,length(x));

uend(1,sites(N+1)) = 1.1*uc; 
%uend(1,:) = d*H(x,tR)/tR;


%a(1,:) = heav(uend(1,sites)-uc);
a(1,:) = prob(uend,1,sites);

G = 1.0/sqrt(4*pi*Diff*tR)*exp(-tR/td)*exp(-(x.^2)./(4*Diff*tR));
fftG = real(fft(fftshift(G)));


refractory = 1;

for p = 2:p_max
    
   refractory = p-1;
   if refractory > refractory_max
      refractory = refractory_max;
   end
    
   
   fftuend = (fft(fftshift(uend(p-1,:))));
   uend(p,:) = uend(p,:) + ifftshift(real(ifft(fftuend.*fftG)))*(step);
  
    
    Hbasic=zeros(size(x));
        for i=1:2*N+1
            if a(p-1,i)==1
                Hbasic = Hbasic + d*H(x-x(sites(i)),tR)/tR;
            end
        end

    
     uend(p,:) = uend(p,:) + Hbasic;
                
      
      %a(p,:) = heav(uend(p,sites)-uc);
      a(p,:)=prob(uend,p,sites);      
      for m=1:refractory
          a(p,:)=a(p,:).*(1-a(p-m,:));
      end
      
     
 end

pcolor(uend)
shading interp
caxis([0 1.2])






