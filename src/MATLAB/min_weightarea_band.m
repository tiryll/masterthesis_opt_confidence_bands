function [clb,coverage] = min_weightarea_band(Omega,Jacobian,weights,nsims,cred,lB)
    
    [Omega_V,Omega_D] = eig(Omega);
    
    s = nsims^(-1/2);
    
    Omega_sqrt = sqrt(Omega_D)*Omega_V';
    Z = randn(nsims,size(Omega,1))*real(Omega_sqrt);

    G=Jacobian;

    M=size(Jacobian,1);

    Vhat=Z*G';

    options = optimoptions('fmincon','Display','notify', 'Algorithm', 'interior-point', ...
        'MaxFunctionEvaluations', 5000, 'MaxIterations', 3000,'SpecifyConstraintGradient', ...
        true,'SpecifyObjectiveGradient', true);

    clb = fmincon(@(param)obj1(param,M,weights),lB,[],[],[],[], [], [], @(param)constr(param,M,Vhat,cred,s,nsims), options);
    
    [c,~,~,~]=constr(clb,M,Vhat,cred,s, nsims);

    coverage=c(1,1);

end

function f = kernf(u)

f = (1+exp(-u)).^(-1);

end


function f = dkernf(u)

y=abs(u);
f = exp(-y)./(1+exp(-y)).^2;

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
    dc =[dc1,dc2];
    
    ceq = [];
    dceq = [];
end
