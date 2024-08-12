
function results = sim_svar_1_vec(alpha,T,p,H,mc_runs)

K=2;

A = [alpha, zeros(1, K*p - 1); 0.5.*repelem((1:p).^(-2), K)];

comp_form = [A;[eye(K*(p-1)), zeros(K*(p-1), K)]];

P = [1, 0; 0.3, sqrt(1 - 0.3^2)];
Sigma_U = P * P';

W = randn(T+p,K,mc_runs);
y = zeros(T + K*p,K,mc_runs);

for t = 1:(T + p)
    y(t + p, :,:) = pagemtimes(A,reshape(permute(y((t + p - 1):-1:t, :,:), ...
        [2,1,3]), [], 1,mc_runs)) + pagemtimes(P,permute(W(t, :,:),[2,1,3]));
end

Phi = zeros(K, K, H + 1);
Theta = zeros(K, K, H + 1);

J = [eye(2), zeros(2, 2*(p - 1))];

for h = 0:H
    Phi(:, :, h+1) = (J * (comp_form^h) * J');
    Theta(:, :, h+1) = Phi(:, :, h + 1) * P;
end

results.y = y((2*p+1):end, :, :);
results.param = [alpha,T,p,H];
results.Theta = Theta;
results.Phi = Phi;

end