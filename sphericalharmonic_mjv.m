%Probability density using legendre polynomials
LL=0;% This is the input L value for the legendre() function
Nthet=200; %parameter to define the number of points at which the 
% polynomials will be evaluated
STEPTHET=pi/Nthet;%Define step size in terms of Nthet
Nphi = 200;
STEPPHI = 2*pi/Nphi;
thet=linspace(0+STEPTHET/2, pi-STEPTHET/2,Nthet);% linspace generates a vector
% of theta values. We cut off the linspace function before pi to avoid
% double counting
phi=linspace(0+STEPPHI/2, 2*pi-STEPPHI/2,Nphi);
[THET,PHI]=meshgrid(thet,phi);%use meshgrid to define input matrices
costhet=cos(THET);% generate input of cosine(theta)
N=legendre(LL,costhet,'norm');% use legendre() to generate a set of
%polynomials dependent on the input LL

figure(1)
r = N(:,:,1)' * N(:,:,1); % r is a matrix of output values

[x,y,z]= sph2cart(THET,PHI,r); %convert from spherical coord to cartesian
%plot with surf
surf(x,y,z);

%Here we generate a the legendre polynomials again with LL = 1
LL=1;

N=legendre(LL,costhet,'norm');%
%We have to change the order of the output matrix to allow 
%easy access to its values. We use permute to do that. 
Nt = permute(N, [2 3 1]);

%Here another figure is plotted
figure(2)
%Pull values from the first 200x200 array that Nt contains
r = Nt(:,:,1)' * Nt(:,:,1); 

%The values of z and -z are plotted, similar to above. 
[x,y,z]= sph2cart(THET,PHI,r); 

surf(x,y,z);
hold on 

surf(x,y,-z)

%%%% Here, I'm just testing out how my code performs for 
%%%% higher order polynomials

% LL=2;
% 
% N=legendre(LL,costhet,'norm');%
% %We have to change the order of the output matrix to allow 
% %easy access to its values. We use permute to do that. 
% Nt = permute(N, [2 3 1]);
% 
% %Here another figure is plotted
% figure(3)
% %Pull values from the first 200x200 array that Nt contains
% r = Nt(:,:,1)' * Nt(:,:,1); 
% 
% %The values of z and -z are plotted, similar to above. 
% [x,y,z]= sph2cart(THET,PHI,r); 
% 
% surf(x,y,z);
% hold on 
% 
% % surf(x,y,-z)
% % r = Nt(:,:,2)' * Nt(:,:,1);
% % surf(x,y,z);
% % 
% % surf(x,y,-z)