
function results = cov_var_NID_vec(Sigma_U, Gamma, d, p)

Sigma_A = kron_vec(pageinv(Gamma), Sigma_U);

K = size(Sigma_U, 1);

D = D_matrix(K);
D_plus = (D' * D) \ D';

Sigma_sigma = 2.*pagemtimes(D_plus,pagemtimes(kron_vec(Sigma_U, Sigma_U),D_plus'));

len_vecA = K * (d + K * p);
len_vechSigma = (1/2) * K * (K + 1);
len_Omega = len_vecA + len_vechSigma;

Omega = zeros(len_Omega, len_Omega,size(Gamma,3));

Omega(1:len_vecA, 1:len_vecA,:) = Sigma_A.*ones(len_vecA,len_vecA,size(Gamma,3));
Omega((len_vecA+1):len_Omega, (len_vecA+1):len_Omega,:) = Sigma_sigma.*ones(len_vechSigma,len_vechSigma,size(Gamma,3));

results.Omega = Omega;
results.Sigma_A = Sigma_A;
results.Sigma_U = Sigma_U;
results.const = d;

end

function out = kron_vec(A,B)

    out =zeros(size(A,1)*size(B,1),size(A,2)*size(B,2),size(A,3));

    for r=1:size(A,1)
        for c=1:size(A,2)
            
            out((r-1)*size(B,1)+1:r*size(B,1), ...
                (c-1)*size(B,2)+1:c*size(B,2),:)=A(r,c,:).*B;
        end
    end

end


function D = D_matrix(K)
  
    auxmat1 = reshape(1:K^2,K,K)-[0, (2:K).*((2:K)-1)/2];
    auxmat2= auxmat1 - tril(auxmat1')' + tril(auxmat1)';
    
    vec=reshape(auxmat2,[],1);
    
    D=1*(vec==1:K*(K+1)/2);
end