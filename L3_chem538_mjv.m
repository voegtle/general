%Using the legendre() function to generate a set of associated legendre
%polynomials
Nthet=100; %parameter to define the number of points at which the 
% polynomials will be evaluated
STEPLG=pi/Nthet;%Define step size in terms of Nthet
LL=1; %Specifies the order of the polynomial
thet=linspace(0+STEPLG/2, pi-STEPLG/2,Nthet);% linspace generates a vector
% of theta values. We cut off the linspace function before pi to avoid
% double counting
costhet=cos(thet);
sinthet=sin(thet);
LLthet=legendre(LL,costhet,'norm');% generate associated leg. polynomials

figure(1);
polarplot(thet,LLthet'*LLthet/(2*pi));% we're using polar plot in this case
ax=gca;
ax.ThetaZeroLocation='top';
ax.ThetaDir='clockwise';

OVLG=LLthet.*sinthet*LLthet'*STEPLG;% test for orthogonality with overlap



%Associated Laguerre Functions
STEP3=0.05;
x3=0:STEP3:50;
n=5;L=2; a0=1;
LAG=laguerreL(n,L,x3); %Use laguerreL() to define polynomials over 
%specified x range 
figure (2);
plot(x3,LAG);%Plot of Laguerre Polynomials v. X
ylim([-5,10]);

Hfn= @(n,L,x,a0)(2*x/(n*a0)).^L.*laguerreL(n-1-L,2.*L+1,(2*x/(n*a0)))...
    .*exp(-2*x/(n.*a0)/2).*sqrt((2/(n.*a0)).^3.*...
    factorial(n-L-1)/(2.*n.*factorial(n+L)));
% Definition of a custom function. Generates Lag. Polynomials with an 
% appropriate weighting

wfn=Hfn(n,L,x3,a0); %uses Hfn to define wfn
n=5;L=2;

figure(3);
for jj=1+L:n;
% this loop iterates over n values from 1 + L to n, then plots the result
    wfn=Hfn(jj,L,x3,a0);
    
    plot(x3,wfn);
    legstr='L+1';
    legend(legstr)
    xlim([0 50]);
    hold on;
end
legend('n=3','n=4','n=5');
hold off

%Normalize
Norm=x3.^2.*Hfn(4,3,x3,a0)*Hfn(4,3,x3,a0)'*STEP3;

Hfn(5,3,100,a0)

x3=0:.1:100;
%Overlap
Overlap=x3.^2.*Hfn(4,3,x3,a0)*Hfn(5,3,x3,a0)'*STEP3;