function results = sim_svar_1_irf_coeff(rho,p,H)

d=0;
K=2;

len_vecA = K * (d + K * p);
len_vechSigma = (1/2) * K * (K + 1);
len_Omega = len_vecA + len_vechSigma;

A_true = [rho, zeros(1, 2*p - 1); 0.5.*repelem((1:p).^(-2), 2)];
A_true_comp_form = [A_true; [eye(2*(p-1)), zeros(2*(p-1), 2)]];

P = [1, 0; 0.3, sqrt(1 - 0.3^2)];
Sigma_U_true = P * P';

Phi_true = zeros(2, 2, H + 1);
Theta_true = zeros(2, 2, H + 1);

Phi_true(:, :, 1) = eye(2);
Theta_true(:, :, 1) = P;

J = [eye(2), zeros(2, 2*(p - 1))];

for h = 1:H
    Phi_true(:, :, h + 1) = J * (A_true_comp_form^h) * J';
    Theta_true(:, :, h + 1) = Phi_true(:, :, h + 1) * P;
end

% true Omega and Sigma

if abs(rho)<1

sigma_U_true_comp_form_vec=reshape([Sigma_U_true,zeros(K,K*(p-1)); ...
    zeros(K*(p-1),K*p)],[],1);

Rho_0_Y_comp_form_vec=(eye((K*p)^2)-kron(A_true_comp_form,A_true_comp_form))\sigma_U_true_comp_form_vec;
Rho_0_Y_comp_form=reshape(Rho_0_Y_comp_form_vec,[K*p,K*p]);

NID_covs=cov_var_NID(Sigma_U_true,Rho_0_Y_comp_form,d,p);

Omega_true=NID_covs.Omega;
Sigma_A_true=NID_covs.Sigma_A;

irf_delta_cov=delta_irf_chol(Phi_true,Omega_true,Sigma_A_true,A_true_comp_form,P,d,H);

Sigma_true=irf_delta_cov.Sigma_Theta;
Jacobian_true=irf_delta_cov.Jacobian;

else
    Omega_true=NaN(len_Omega,len_Omega);
    Sigma_true=NaN(K^2*(H+1),K^2*(H+1));
    Jacobian_true=NaN(K^2*(H+1),len_Omega);
end

results.A_true_comp_form = A_true_comp_form;
results.Sigma_U_true = Sigma_U_true;
results.Theta_true = Theta_true;
results.Phi_true = Phi_true;
results.Omega_true=Omega_true;
results.Sigma_true=Sigma_true;
results.Jacobian_true=Jacobian_true;

end
