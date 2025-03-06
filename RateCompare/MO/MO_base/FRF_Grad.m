function egrad = FRF_Grad(x,Heff_l, phi_l, W_l, T_l, Nt, Nt_rf,K)


x = reshape(x,Nt,Nt_rf);
for k = 1:K
    M_l(:,:,k) = T_l(:,:,k)^(-1)+phi_l(1,k)^(-1)*...
        Heff_l(:,:,k)*x*(x'*x)^(-1)*x'*Heff_l(:,:,k)';
    MM_l(:,:,k) = phi_l(1,k)^(-1)*Heff_l(:,:,k)'*...
        M_l(:,:,k)^(-2)*Heff_l(:,:,k)*x;
end
M_sum = sum(MM_l,3);
egrad = (x*(x'*x)^(-1)*x'-eye(Nt))*M_sum*(x'*x)^(-1);
egrad = egrad(:);