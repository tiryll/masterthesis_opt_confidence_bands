
function results = cov_var_NID(Sigma_U, Gamma, d, p)

Sigma_A = kron(inv(Gamma), Sigma_U);

K = size(Sigma_U, 1);

D = D_matrix(K);
D_plus = (D' * D) \ D';

Sigma_sigma = 2 * D_plus * kron(Sigma_U, Sigma_U) * D_plus';

len_vecA = K * (d + K * p);
len_vechSigma = (1/2) * K * (K + 1);
len_Omega = len_vecA + len_vechSigma;

Omega = zeros(len_Omega, len_Omega);

Omega(1:len_vecA, 1:len_vecA) = Sigma_A;
Omega((len_vecA+1):len_Omega, (len_vecA+1):len_Omega) = Sigma_sigma;

results.Omega = Omega;
results.Sigma_A = Sigma_A;
results.Sigma_U = Sigma_U;
results.const = d;

end


function D = D_matrix(K)
  
    auxmat1 = reshape(1:K^2,K,K)-[0, (2:K).*((2:K)-1)/2];
    auxmat2= auxmat1 - tril(auxmat1')' + tril(auxmat1)';
    
    vec=reshape(auxmat2,[],1);
    
    D=1*(vec==1:K*(K+1)/2);
end