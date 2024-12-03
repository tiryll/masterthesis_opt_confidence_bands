function results=cov_var_HAC(Y,U,Sigma_U,Gamma,d,p)

K = size(Sigma_U, 1);

D = D_matrix(K);
D_plus = (D' * D) \ D';

Sigma_A=kron(inv(Gamma(1+d:end,1+d:end)),Sigma_U);
Sigma_sigma = D_plus * kron(Sigma_U, Sigma_U) * D_plus';

J_inv=[0.5*Sigma_A,zeros(size(Sigma_A,1),size(Sigma_sigma,2));
    zeros(size(Sigma_sigma,1),size(Sigma_A,2)),Sigma_sigma];

X = Y((p+1):end, :);
T = size(X, 1);
Z = zeros(T, p*K);

for ii = 1:p
    Z(:, ((ii-1)*K+1):(ii*K)) = Y((p+1-ii):(end-ii), :);
end

Z=Z';

len_Ups=size(Sigma_A,1)+size(Sigma_sigma,1);

Ups=NaN(len_Ups,1,T-1);

Ups(1:size(Sigma_A,1),:,:)=pagemtimes(2*reshape(kron(Z(:,1:end-1),eye(K)),[p*K^2,K,T-1]), ...
    reshape(Sigma_U\U(:,2:end),[K,1,T-1]));

for t=1:T-1
    
    Ups(size(Sigma_A,1)+1:len_Ups,:,t)=D_matrix(K)'*(reshape(inv(Sigma_U),[K^2,1]) - reshape(Sigma_U\(U(:,t+1)*U(:,t+1)')*inv(Sigma_U),[K^2,1]));

end

Ups=reshape(Ups,[len_Ups,T-1]);

aux_var_Ups=fit_rvar_lagselect(Ups',0,1,10); % selects auxiliary VAR model by AIC criterion with p_max=10

A_aux=aux_var_Ups.A;
Sigma_U_aux=aux_var_Ups.Sigma_U_MLE;
p_aux=aux_var_Ups.p;

char_pol_aux=sum(cat(3,eye(len_Ups),reshape(A_aux,len_Ups,len_Ups,p_aux)),3);

B=char_pol_aux \ Sigma_U_aux / char_pol_aux';
results.V=J_inv*B*J_inv';
results.Sigma_A=Sigma_A;

end

function D = D_matrix(K)
  
    auxmat1 = reshape(1:K^2,K,K)-[0, (2:K).*((2:K)-1)/2];
    auxmat2= auxmat1 - tril(auxmat1')' + tril(auxmat1)';
    
    vec=reshape(auxmat2,[],1);
    
    D=1*(vec==1:K*(K+1)/2);
end