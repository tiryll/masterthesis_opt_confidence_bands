
function supt_quant = supt_irf_delta(Sigma_irf,index,cred,nsims)

Sigma_select = Sigma_irf(index,index);

nzero_entries=diag(Sigma_select)>0;
Sigma_select_nzero = Sigma_select(nzero_entries,nzero_entries);

stdev_Sigma = sqrt(diag(Sigma_select_nzero));

Rho=(stdev_Sigma.^(-1)).*Sigma_select_nzero.*(stdev_Sigma.^(-1))';

[Rho_E,Rho_D] = eig(Rho);

sim_supt=max(abs(Rho_E*sqrt(Rho_D)*randn(sum(nzero_entries),nsims)));

supt_quant=quantile(sim_supt,cred);

end