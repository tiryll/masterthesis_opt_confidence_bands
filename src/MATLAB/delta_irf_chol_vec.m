function results = delta_irf_chol_vec(Phi,Omega, Sigma_A,A_comp_form,P,d,H)

A = A_comp_form;

K = size(P, 1);
p = size(A, 2) / K;
mc_runs=size(Phi,4);

len_vecirf = K^2;
len_vecA = p * K^2;
len_vechSigma = K * (K + 1) / 2;

if d >= 1
    Omega_irf = Omega(K * d + 1:end, K * d + 1:end,:);
    Sigma_A_irf = Sigma_A(K * d + 1:end, K * d + 1:end,:);
else
    Omega_irf = Omega;
    Sigma_A_irf = Sigma_A;
end

L_mat = L_matrix(K);
K_mat = K_matrix(K);

Sigma_Phi=zeros(len_vecirf*H,len_vecirf*H,mc_runs);
Sigma_Theta=zeros(len_vecirf*(H+1),len_vecirf*(H+1),mc_runs);
Jacobian_Theta=zeros(len_vecirf*(H+1),len_vecA+len_vechSigma,mc_runs);

for r=1:mc_runs

    Jacobian_Phi_1_H_r = zeros(len_vecirf, len_vecA, H);

    for h = 1:H
        for m = 0:(h-1)
            A_pow = (A(:,:,r)')^(h - 1 - m);
            Jacobian_Phi_1_H_r(:, :, h) = Jacobian_Phi_1_H_r(:, :, h) + kron(A_pow(1:K, :), Phi(:, :, m + 1,r));
        end
    end

    Jacobian_Phi_r = reshape(permute(Jacobian_Phi_1_H_r, [2, 1, 3]), len_vecA, len_vecirf*H)';
    Sigma_Phi(:,:,r) = Jacobian_Phi_r * Sigma_A_irf(:,:,r) * Jacobian_Phi_r';

    Jacobian_Theta_dA_r = cat(3,zeros(len_vecirf,len_vecA), ...
    pagemtimes(kron(P(:,:,r)',eye(K)).*ones(len_vecirf,len_vecirf,H),Jacobian_Phi_1_H_r));

    Jacobian_Theta_dSigma_r = zeros(len_vecirf,len_vechSigma,H+1);

    H_mat = L_mat' / (L_mat * (eye(K^2) + K_mat) * kron(P(:,:,r), eye(K)) * L_mat');

    for h=1:H+1
        Jacobian_Theta_dSigma_r(:,:,h) = kron(eye(K),Phi(:,:,h,r))*H_mat;
    end

    Jacobian_Theta_r = reshape(permute(cat(2,Jacobian_Theta_dA_r,Jacobian_Theta_dSigma_r),[2 1 3]), ...
    len_vecA+len_vechSigma,len_vecirf*(H+1))';

    Sigma_Theta(:,:,r) = Jacobian_Theta_r * Omega_irf(:,:,r) * Jacobian_Theta_r';
    Jacobian_Theta(:,:,r)=Jacobian_Theta_r;

end

results.Sigma_Phi = Sigma_Phi;
results.Sigma_Theta = Sigma_Theta;
results.nvars = K;
results.horizon = H;
results.Jacobian = Jacobian_Theta;
end

function L_mat =L_matrix(K)

    auxmat1 = reshape(1:K^2,K,K)-[0, (2:K).*((2:K)-1)/2];
    auxmat2= auxmat1 - tril(auxmat1',-1)';
    
    vec0=reshape(auxmat2,1,[]);
    
    L_mat=1*((1:K*(K+1)/2)'==vec0);


end

function K_mat=K_matrix(K)

    auxmat1=reshape(1:K^2,K,K);
    auxmat2=auxmat1';

    K_mat=1*(reshape(auxmat1,[],1)==reshape(auxmat2,1,[]));

end