
function xmin = quantiles_wmeanLoss_sym(Sigma_Theta,index,weights,nsims,cred,lB)
    
    options = optimoptions('fmincon','Display','notify', 'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 50000, 'MaxIterations', 30000,'SpecifyConstraintGradient',true,'SpecifyObjectiveGradient', true);
    
    Cov=Sigma_Theta(index,index);

    nzero_entries = diag(Cov)>0;
    K=sum(nzero_entries);

    Cov_nzero=Cov(nzero_entries,nzero_entries);

    [Cov_V,Cov_D] = eig(Cov_nzero);
    Cov_D(Cov_D<0)=0;

    Cov_sqrt = sqrt(real(Cov_D))*real(Cov_V)';

    hn = nsims^(-1/2);
    
    Zs = randn(nsims,K)*Cov_sqrt;

    xmin = fmincon(@(param)obj1(param,K,weights), lB, [], [], [], [], [], [], @(param)constr1(param, Zs,K,cred,hn,nsims), options);

end

function f = kernf(u)

f = ((3/4)*( (u+1) - (1/3)*(u.^3+1) )).*(u>=-1).*(u<=1) + (u>1);

end


function f = dkernf(u)

f = ((3/4)*( 1 - u.^2 )).*(u>=-1).*(u<=1);

end

function [f,df] = obj1(param,K,weights)

    lB = param(1:K,1);
    uB = -lB;


    f = mean( (uB - lB).*weights);

    df =  -2.*weights/K;
end

function [c, ceq, dc, dceq] = constr1(param, Zs, K,cred,hn, ndraws)

    lB = param(1:K,1);
    uB = -lB;
   
    arg1 = Zs-ones(ndraws,1)*(lB');
    arg2 = ones(ndraws,1)*(uB')- Zs;
    
    cmat = [kernf(arg1/hn) kernf(arg2/hn)];
    
    dcmat = [-dkernf(arg1/hn)/hn -dkernf(arg2/hn)/hn];
     
    indall = (1:2*K)'; 
    
    dc1 = zeros(K,1);
    for j = 1:K
        
    temp1 = (indall~=j);
    temp2 = (indall~=(j+K));
    
    dc1(j,1) = -( mean(prod(cmat(:,indall(temp1)),2).*dcmat(:,j) + prod(cmat(:,indall(temp2)),2).*dcmat(:,j+K)));
    
    end
    
    c1 = cred - mean((prod(kernf(arg1/hn),2).*prod(kernf(arg2/hn),2)));
       
    c2 = 0.1 - (uB - lB);
    dc2 = 2*eye(K);

	c = [c1;c2];
    dc = [dc1 dc2];
    
    ceq = [];
    dceq = [];
end