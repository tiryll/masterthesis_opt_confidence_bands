
function results=var_arch_test(U,q,alpha)

    K=size(U,1);
    T=size(U,2);

    len_vech_outer=1/2*K*(K+1);
    vech_outer_U=NaN(len_vech_outer,T);

    for t=1:T

        outer_U=U(:,t)*U(:,t)';
        vech_outer_U(:,t)=L_matrix(K)*reshape(outer_U,[],1);

    end

    Sigma_0=1/T*(vech_outer_U*vech_outer_U');

    aux_var=fit_var(vech_outer_U',q,1);
    Sigma_vech=aux_var.Sigma_U_MLE;

    MARCH_LM_stat=1/2*T*K*(K+1)-T*trace(Sigma_vech/Sigma_0);

    df=q*((K^2+K)^2)/4;

    results.p_val=chi2cdf(MARCH_LM_stat,df);
    results.crit_val=chi2inv(1-alpha,df);
    results.desicion=(results.crit_val<MARCH_LM_stat);
    results.test_stat=MARCH_LM_stat;

end

function L =L_matrix(K)

    auxmat1 = reshape(1:K^2,K,K)-[0, (2:K).*((2:K)-1)/2];
    auxmat2= auxmat1 - tril(auxmat1',-1)';
    
    vec0=reshape(auxmat2,1,[]);
    
    L=1*((1:K*(K+1)/2)'==vec0);


end