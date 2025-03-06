function cost = FRF_Cost(x,Heff_l, phi_l, W_l, T_l, Nt, Nt_rf,K)



x = reshape(x,Nt,Nt_rf);

for k = 1:K
    cost(k) = trace((T_l(:,:,k)^(-1)+phi_l(1,k)^(-1)*...
        Heff_l(:,:,k)*x*(x'*x)^(-1)*x'*Heff_l(:,:,k)')^(-1));
end
cost = sum(cost);

end