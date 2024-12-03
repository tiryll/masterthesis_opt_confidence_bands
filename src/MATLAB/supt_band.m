
function supt_quant = supt_band(Sigma,cred,Ns)

M=size(Sigma,1);

stdev_Sigma = sqrt(diag(Sigma));

Rho=(stdev_Sigma.^(-1)).*Sigma.*(stdev_Sigma.^(-1))';

[Rho_E,Rho_D] = eig(Rho);

sim_supt=max(abs(Rho_E*sqrt(Rho_D)*randn(M,Ns)));

supt_quant=quantile(sim_supt,cred);

end