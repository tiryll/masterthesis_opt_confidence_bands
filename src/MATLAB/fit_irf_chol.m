
function results = fit_irf_chol(A,Sigma_U,U,p,d,H)

K = size(U, 1);

if d == 1
    A_irf = A(:, 2:end);
else
    A_irf = A;
end

A_comp_form = [A_irf; [eye(K*(p-1)), zeros(K*(p-1), K)]];

P = chol(Sigma_U, 'lower');
W = P \ U;

J = [eye(K), zeros(K, K*(p-1))];

Phi = zeros(K, K, H+1);
Theta = zeros(K, K, H+1);

for h = 0:H
    Phi(:, :, h+1) = J * (A_comp_form^h) * J';
    Theta(:, :, h+1) = Phi(:, :, h+1) * P;
end

results.Theta = Theta;
results.Phi = Phi;
results.A_comp_form = A_comp_form;
results.W = W;
results.P = P;
results.J = J;
results.H = H;

end
