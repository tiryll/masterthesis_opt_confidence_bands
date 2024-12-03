
function results = fit_irf_chol_vec(A,Sigma_U,U,p,d,H)

K = size(U, 1);
mc_runs=size(A,3);
if d == 1
    A_irf = A(:, 2:end,:);
else
    A_irf = A;
end

A_comp_form = cat(1,A_irf,[eye(K*(p-1)), zeros(K*(p-1), K)].*ones(K*(p-1),K*p,mc_runs));

P=zeros(K,K,mc_runs);

J = [eye(K), zeros(K, K*(p-1))];

Phi = zeros(K, K, H+1,mc_runs);
Theta = zeros(K, K, H+1,mc_runs);

for r=1:mc_runs
P(:,:,r) = chol(Sigma_U(:,:,r), 'lower');
    for h = 0:H
        Phi(:, :, h+1,r) = J * (A_comp_form(:,:,r)^h) * J';
        Theta(:, :, h+1,r) = Phi(:, :, h+1,r)*P(:,:,r);
    end
end

W = pagemtimes(pageinv(P),U);

results.Theta = Theta;
results.Phi = Phi;
results.A_comp_form = A_comp_form;
results.W = W;
results.P = P;
results.J = J;
results.H = H;

end
