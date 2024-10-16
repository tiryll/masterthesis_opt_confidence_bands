function clb = quantiles_wmeanLoss_sym3(Omega,Jacobian,index,weights,nsims,cred,lB)
    
    options = optimoptions('fmincon','Display','notify', 'Algorithm', 'interior-point', ...
        'MaxFunctionEvaluations', 50000, 'MaxIterations', 30000,'SpecifyConstraintGradient', ...
        true,'SpecifyObjectiveGradient', true);

    [Omega_V,Omega_D] = eig(Omega);
    
    s = nsims^(-1/2);
    
    Omega_sqrt = sqrt(Omega_D)*Omega_V';
    Z = randn(nsims,size(Omega,1))*Omega_sqrt;

    G=Jacobian(index,:);

    M=length(index);

    Vhat=Z*G';

    clb_opt = fmincon(@(param)obj1(param,M,weights),lB,[],[],[],[], [], [], @(param)constr(param,M,Vhat,cred,s,nsims), options);
    
    % ensure validity 
    G_scaled = G.*((clb_opt(1:M,1).^(-1))*ones(1,size(Omega,1)));
    b = quantile(max(abs(G_scaled*Z'))',cred);

    clb = clb_opt(1:M,1).*b;

end

function f = kernf(u)

f = (1+exp(-u)).^(-1);

end


function f = dkernf(u)

f = exp(-u-2.*log(1+exp(-u)));

end

function f=ddkernf(u)

f=exp(1+exp(-u)).^(-3)+exp(2-u-3.*log(1+exp(-u)))+exp(4-2.*u-3.*log(1+exp(-u)))+exp(-3.*u-3.*log(1+exp(-u)))-1;

end

function [f,df]= obj1(param,M,weights)

    lB = param(1:M,1);
    uB = -lB;


    f = sum( (uB - lB).*weights);

    df =  -2.*weights;
end

function [c, ceq, dc, dceq] = constr(param,M,Vhat,cred,s, nsims)
   
    lB = param(1:M,1);
    uB = -lB;

    arg1 = Vhat-ones(nsims,1)*(lB');
    arg2 = ones(nsims,1)*(uB')- Vhat;
    
    cmat = [kernf(arg1/s) kernf(arg2/s)];
    
    dcmat = [-dkernf(arg1/s)/s -dkernf(arg2/s)/s];
     
    indall = (1:2*M)'; 
    
    dc1 = zeros(M,1);
    for j = 1:M
        
    temp1 = (indall~=j);
    temp2 = (indall~=(j+M));
    
    dc1(j,1) = -( mean(prod(cmat(:,indall(temp1)),2).*dcmat(:,j) + prod(cmat(:,indall(temp2)),2).*dcmat(:,j+M)));
    
    end
    
    c1 = cred - mean((prod(kernf(arg1/s),2).*prod(kernf(arg2/s),2))); % coverage constraint
       
    c2 = 0.1 - (uB - lB); %non-negativity constraint
    dc2 = 2*eye(M);

	c = [c1;c2];
    dc = [dc1 dc2];
    
    ceq = [];
    dceq = [];
end

function h = hessian(param,M,Vhat,cred,s, nsims)

    lB = param(1:M,1);
    uB = -lB;

    arg1 = Vhat-ones(nsims,1)*(lB');
    arg2 = ones(nsims,1)*(uB')- Vhat;
    
    cmat = [kernf(arg1/s) kernf(arg2/s)];
    
    dcmat = [-dkernf(arg1/s)/s -dkernf(arg2/s)/s]; % matrix of first derivatives of F
    ddcmat = [ddkernf(arg1/s)/s^2 ddkernf(arg2/s)/s^2]; % matrix of second derivatives of F
    
    indall = (1:2*M)'; 
    
    ddc1 = zeros(M,M);

    for j = 1:M

        temp1_j = (indall~=j);
        temp2_j = (indall~=(j+M));
        
        for m=1:M

        temp1_m = (indall~=m);
        temp2_m = (indall~=(m+M));

        temp1=temp1_j & temp1_m;
        temp2=temp2_j & temp2_m;
    
            if m~=j
                ddc1(j,1) = -( mean(prod(cmat(:,indall(temp1)),2).*dcmat(:,j).*dcmat(:,) + prod(cmat(:,indall(temp2)),2).*dcmat(:,j+M)));
        
        end
    end
    
    c1 = cred - mean((prod(kernf(arg1/s),2).*prod(kernf(arg2/s),2))); % coverage constraint
  
end