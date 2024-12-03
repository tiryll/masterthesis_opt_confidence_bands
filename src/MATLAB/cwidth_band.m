
function cwidth_quant = cwidth_band(Sigma,cred,Ns)

    M=size(Sigma,1);

    [Sigma_V,Sigma_D] = eig(Sigma);
    
    Sigma_D(Sigma_D<0)=0;
    Sigma_sqrt = sqrt(real(Sigma_D))*real(Sigma_V)';
    
    Zs = randn(Ns,M)*Sigma_sqrt;

    sim_cwidth=max(abs(Zs),[],2);
    cwidth_quant=quantile(sim_cwidth,cred);

end