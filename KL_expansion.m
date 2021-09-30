% 3D KL-expansion simulation %
% The parameter setting and covariance matrix construction are the same as code above (Matrix decomposition method, Subsection 2.4)

% Here we apply the SVD method,
% K-L transform is a common orthogonal transform used for data compression and dimensionality reduction
%Orthogonal transformation means that when we analyze the components of one frequency, by orthogonal transformation,
%the components of the analysis will not be mixed with the components of other frequencies
%(for example, if a wave is a combination of low and high frequencies, we want to decompose the two waveforms by transformation).
%How do you do an orthogonal transformation, suppose you have an orthogonal basis, and then you use that orthogonal basis to decompose a signal,
%and you have coefficients in front of it, and you convert x into these coefficients
%% SVD,
% Parameter Setting %
nx = 10;
ny = 10;
nz = 10;

Lx = 1;
Ly = 1;
Lz = 1;

dx = Lx / nx;
dy = Ly / ny;
dz = Lz / nz;

covMat=zeros(nx*ny);  %zeros(n) returns an n by n all-zero matrix

% Matern covariance, with parameters 
nu=1;
lambda=1;
sigma2=1;

% Measure distance between ii and jj, construct covariance matrix
for ii=1:nx*ny*nz
    for jj=1:nx*ny*nz
        pageii = floor(ii/(nx*ny)); %Floor rounds down and returns the largest integer less than or equal to a given number
        premi   = rem(ii,(nx*ny));  %The remainder of rem divided by r equals rem of a,b returns the remainder of a divided by b
        rowii  = floor(premi/nx);
        colii  = rem(premi,nx);
        
        pagejj = floor(jj/(nx*ny));
        premj   = rem(jj,(nx*ny));
        rowjj  = floor(premj/nx);
        coljj  = rem(premj,nx);
        
        z=abs((pagejj-pageii))*dz;
        y=abs((rowjj-rowii))*dy;
        x=abs((coljj-colii))*dx;  %ii = nx*ny*pageii + nx*rowii + colii
        d=(z^2+y^2+x^2)^(1/2);
        if (d==0)
            covMat(ii,jj)=sigma2;
        else
         %K=besselk(nu,z) computes the second type of modified BESselk function for each element in array Z
        covMat(ii,jj)=sigma2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*d/lambda)^(nu)*besselk(nu,sqrt(2*nu)*d/lambda); 
        end
    end
end

Nkl = 300;
%[V,D] = eigs(___) returns the diagonal matrix D containing the eigenvalues along the main diagonal, 
%and the matrix V containing the corresponding eigenvectors in the columns of the latter
[eigenVec,eigenVal]=eigs(covMat,Nkl);
%%
eigenVal=diag(eigenVal);  %returns the column vector of the major diagonal element of A
%%
xi = randn(Nkl,1); %Returns an NK1*1 matrix of normally distributed random numbers
K = eigenVec * (sqrt(eigenVal).*xi);

Kr = reshape(K,[nx ny nz]); %Refactoring it to nx times ny times nz
Kre = exp(Kr);
