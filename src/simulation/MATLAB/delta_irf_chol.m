function results = delta_irf_chol(Phi,Omega, Sigma_A,A_comp_form,P,d,H)

A = A_comp_form;

K = size(P, 1);
p = size(A, 2) / K;

len_vecirf = K^2;
len_vecA = p * K^2;
len_vechSigma = K * (K + 1) / 2;

if d >= 1
    Omega_irf = Omega(K * d + 1:end, K * d + 1:end);
    Sigma_A_irf = Sigma_A(K * d + 1:end, K * d + 1:end);
else
    Omega_irf = Omega;
    Sigma_A_irf = Sigma_A;
end

Jacobian_Phi_1_H = zeros(len_vecirf, len_vecA, H);

for h = 1:H
    for m = 0:(h-1)
        A_pow = (A')^(h - 1 - m);
        Jacobian_Phi_1_H(:, :, h) = Jacobian_Phi_1_H(:, :, h) + kron(A_pow(1:K, :), Phi(:, :, m + 1));
    end
end

Jacobian_Phi = reshape(permute(Jacobian_Phi_1_H, [2, 1, 3]), len_vecA, len_vecirf*H)';
Sigma_Phi = Jacobian_Phi * Sigma_A_irf * Jacobian_Phi';

L = L_matrix(K);
K_mat = K_matrix(K);

H_mat = L' / (L * (eye(K^2) + K_mat) * kron(P, eye(K)) * L');

Jacobian_Theta_dA = cat(3,zeros(len_vecirf,len_vecA), ...
    pagemtimes(kron(P',eye(K)).*ones(len_vecirf,len_vecirf,H),Jacobian_Phi_1_H));

Jacobian_Theta_dSigma = zeros(len_vecirf,len_vechSigma,H+1);

for h=1:H+1

    Jacobian_Theta_dSigma(:,:,h) = kron(eye(K),Phi(:,:,h))*H_mat;

end

Jacobian_Theta = reshape(permute(cat(2,Jacobian_Theta_dA,Jacobian_Theta_dSigma),[2 1 3]), ...
    len_vecA+len_vechSigma,len_vecirf*(H+1))';

Sigma_Theta = Jacobian_Theta * Omega_irf * Jacobian_Theta';

results.Sigma_Phi = Sigma_Phi;
results.Sigma_Theta = Sigma_Theta;
results.nvars = K;
results.horizon = H;

end

function L =L_matrix(K)

    auxmat1 = reshape(1:K^2,K,K)-[0, (2:K).*((2:K)-1)/2];
    auxmat2= auxmat1 - tril(auxmat1',-1)';
    
    vec0=reshape(auxmat2,1,[]);
    
    L=1*((1:K*(K+1)/2)'==vec0);


end

function K_mat=K_matrix(K)

    auxmat1=reshape(1:K^2,K,K);
    auxmat2=auxmat1';

    K_mat=1*(reshape(auxmat1,[],1)==reshape(auxmat2,1,[]));

end