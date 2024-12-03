
function results=fit_rvar_lagselect(X,d,p_min,p_max)

T=size(X,1);
K=size(X,2);

p_opt=(p_min:p_max)';

AIC_temp=NaN(size(p_opt,1),1);

for p=p_min:p_max

    var_mod_temp=fit_var(X,p,d);
    Sigma_U=var_mod_temp.Sigma_U_MLE;

    AIC_temp(p,1)=log(det(Sigma_U))+(K+p*K^2)*2/T;

end

[~,pstar_ind]=min(AIC_temp);

results=fit_var(X,p_opt(pstar_ind),d);

end