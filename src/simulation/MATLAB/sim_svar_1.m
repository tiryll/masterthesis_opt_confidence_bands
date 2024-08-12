
function results = sim_svar_1(alpha,T,p,H)

A = [alpha, zeros(1, 2*p - 1); 0.5.*repelem((1:p).^(-2), 2)];
comp_form = [A; [eye(2*(p-1)), zeros(2*(p-1), 2)]];

P = [1, 0; 0.3, sqrt(1 - 0.3^2)];
Sigma_U = P * P';

W = randn(T+p,2);
y = zeros(T + 2*p, 2);

for t = 1:(T + p)
    y(t + p, :) = A * reshape(y((t + p - 1):-1:t, :)', [], 1) + P * W(t, :)';
end

Phi = zeros(2, 2, H + 1);
Theta = zeros(2, 2, H + 1);

Phi(:, :, 1) = eye(2);
Theta(:, :, 1) = P;

J = [eye(2), zeros(2, 2*(p - 1))];

for h = 1:H
    Phi(:, :, h + 1) = J * (comp_form^h) * J';
    Theta(:, :, h + 1) = Phi(:, :, h + 1) * P;
end

results.y = y((2*p+1):end, :);
results.coefficients = A;
results.sigma = Sigma_U;
results.structural_IRF = Theta;
results.reducedform_IRF = Phi;

end