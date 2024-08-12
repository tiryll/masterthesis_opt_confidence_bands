clc

clear  

options = optimset('display','iter', 'algorithm', 'interior-point', 'MaxFunEvals', 50000, 'MaxIter', 30000, 'Display', 'iter','GradConstr','on','GradObj', 'on' );

rng(4013)

msim = 1; % Note: We ran the simulations on the UW cluster and ran each simulations separately using a different seed.

% Grid of x for finding optimal band
ng = 30;
xgo = (0:ng)'/ng;
 
% Sample size  
n = 500;

% Coverage rate
cr = 0.90; 
 
% Number of simulations for sup-t
nsims = 50000;

% True beta
beta = [1;0.25; 0.2; 0; -0.1]; 
 
sx = size(xgo,1);
xmat = [ones(size(xgo,1),1) sqrt(3)*(2*xgo-1) sqrt(5)*0.5*(3*(2*xgo-1).^2-1) sqrt(7)*0.5*(5*(2*xgo-1).^3-3*(2*xgo-1)) sqrt(9)*1/8*(35*(2*xgo-1).^4-30*(2*xgo-1).^2+3)];

% Input for approximating probabilities (needs normalized vectors)
normx = sqrt(sum(xmat'.^2))';
xmatnorm = xmat'./(ones(size(xmat,2),1)*normx');

ngl = 300;
xgol = (0:ngl)'/ngl;
xmatl = [ones(size(xgol,1),1) sqrt(3)*(2*xgol-1) sqrt(5)*0.5*(3*(2*xgol-1).^2-1) sqrt(7)*0.5*(5*(2*xgol-1).^3-3*(2*xgol-1)) sqrt(9)*1/8*(35*(2*xgol-1).^4-30*(2*xgol-1).^2+3)];
normxl = sqrt(sum(xmatl'.^2))';
xmatnorml = xmatl'./(ones(size(xmatl,2),1)*normxl');

c3w = zeros(ngl,1);
c4w = zeros(ngl,1);

output_mat = zeros(ngl+1,5,msim);
covp = zeros(msim,5);


xtemp = betarnd(2,2,n,1);
xmatr = [ones(size(xtemp,1),1) sqrt(3)*(2*xtemp-1) sqrt(5)*0.5*(3*(2*xtemp-1).^2-1) sqrt(7)*0.5*(5*(2*xtemp-1).^3-3*(2*xtemp-1)) sqrt(9)*1/8*(35*(2*xtemp-1).^4-30*(2*xtemp-1).^2+3)];

sigx = 0.3*(abs(xtemp-0.5)>=0.25) + 3*(abs(xtemp-0.5)<0.25);
    
Y = xmatr*beta + sigx.*randn(n,1);

betah =  (xmatr'*xmatr)\(xmatr'*Y);

% Mean and covariance matrix of the normal dictribution (i.e of (\hat{\beta} - beta))
K = (xmatr'*xmatr)\(xmatr'*diag(( Y - xmatr*betah ).^2)*xmatr)/(xmatr'*xmatr);
  
    
% Sup-t confidence bands using normalized vectors
Z = randn(nsims,size(beta,1))*chol(K);

weights1 = sqrt(diag(xmatnorml'*K*xmatnorml)).^(-1);
xw1 = xmatnorml'.*(weights1*ones(1,size(K,1)));
c1 = quantile(max(abs(xw1*Z'))',cr);
c1w = c1*(weights1).^(-1);

% Adjust (reverse normalization) - Sup-t band is xmat'beta +/- c1w
c1w = c1w.*normxl;

% Minimum area 
    
% Weights for the objective function
oweights = ones(size(xgo,1),1);

ndraws = 50000;
hn = ndraws^(-1/2);
Zs = randn(ndraws,size(beta,1))*chol(K);
    
lB = -2*ones(size(xgo,1),1);
xmin = fmincon(@(param)obj1(param,sx,normx,oweights), lB, [], [], [], [], [], [], @(param)constr1(param, xmatnorm, Zs, sx,cr,hn,ndraws), options);
    
% Adjust (c3 will be close to 0)
weights3 = (xmin(1:size(xgo,1),1)).^(-1);
xw3 = xmatnorm'.*(weights3*ones(1,size(K,1)));
c3 = quantile(max(abs(xw3*Z'))',cr);

c3ws = -xmin(1:size(xgo,1),1).*normx*c3;
    
% Final projection to fill in grid points

for k=1:ngl+1
    [~,c3w(k,1)] = linprog(xmatl(k,:)',[xmat;-xmat],[c3ws;c3ws]);
end
c3w = -c3w;


function [f,df] = obj1(param,sx, normx, oweights)

    lB = param(1:sx,1);
    uB = -lB;


    f = mean( (uB - lB).*normx.*oweights);

    df =  -2.*normx.*oweights/sx;
end
 
 
function [c, ceq, dc, dceq] = constr1(param, xmatnorm, Zs, sx,cr,hn, ndraws)

    lB = param(1:sx,1);
    uB = -lB;
   
    arg1 = (Zs*xmatnorm)-ones(ndraws,1)*(lB');
    arg2 = ones(ndraws,1)*(uB')-(Zs*xmatnorm);
    
    cmat = [kernf(arg1/hn) kernf(arg2/hn)];
    
    dcmat = [-dkernf(arg1/hn)/hn -dkernf(arg2/hn)/hn];
     
    indall = (1:2*sx)'; 
    
    dc1 = zeros(sx,1);
    for j = 1:sx
        
    temp1 = (indall~=j);
    temp2 = (indall~=(j+sx));
    
    dc1(j,1) = -( mean(prod(cmat(:,indall(temp1)),2).*dcmat(:,j) + prod(cmat(:,indall(temp2)),2).*dcmat(:,j+sx)));
    
    end
    
    c1 = cr - mean((prod(kernf(arg1/hn),2).*(prod(kernf(arg2/hn),2))));
       
    c2 = 0.1 - (uB - lB);
    dc2 = 2*eye(sx);

	c = [c1;c2];
    dc = [dc1 dc2];
    
    ceq = [];
    dceq = [];
end
 

function f = kernf(u)

f = ((3/4)*( (u+1) - (1/3)*(u.^3+1) )).*(u>=-1).*(u<=1) + (u>1);

end


function f = dkernf(u)

f = ((3/4)*( 1 - u.^2 )).*(u>=-1).*(u<=1);

end