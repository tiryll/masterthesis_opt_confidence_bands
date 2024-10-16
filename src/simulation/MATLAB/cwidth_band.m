
function cwidth_quant = cwidth_band(Sigma_irf,index,cred,nsims)

    Cov=Sigma_irf(index,index);

    nzero_entries = diag(Cov)>0;
 
    Cov_nzero=Cov(nzero_entries,nzero_entries);

    [Cov_V,Cov_D] = eig(Cov_nzero);
    
    Cov_D(Cov_D<0)=0;
    Cov_sqrt = sqrt(real(Cov_D))*real(Cov_V)';
    
    Zs = randn(nsims,sum(nzero_entries))*Cov_sqrt;

    sim_cwidth=max(abs(Zs),[],2);
    cwidth_quant=quantile(sim_cwidth,cred);

end