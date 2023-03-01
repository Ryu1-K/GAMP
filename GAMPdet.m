function y = GAMPdet(RX,SIM,CH,G)
H = [real(CH.H) -imag(CH.H); imag(CH.H) real(CH.H)];
HH = H.*H;
Y = [real(RX.Y); imag(RX.Y)];

x_hat = zeros(SIM.M*2,SIM.K);
psi_x = CH.Es/2*ones(SIM.M*2,SIM.K);
s_hat = zeros(SIM.N*2,SIM.K);

for iter = 1:SIM.niter
    psi_p = HH*psi_x;
    z_hat = Y - H*x_hat + psi_p.*s_hat;
    nu_y  = psi_p+CH.N0/2;
    s_hat = z_hat ./ nu_y;
    
    psi_r = 1./(HH.'*(1./nu_y));
    
    if iter == 1
        r_hat = (x_hat + psi_r.*(H.'*s_hat));
    else
        r_hat = (SIM.eta)*(x_hat + psi_r.*(H'*s_hat))+(1-SIM.eta)*r_hat;
    end
    
    if SIM.asb ==1
        x_hat_ = G.xc*tanh(SIM.alpha(iter)./G.xc.*(r_hat));
    else
        x_hat_ = G.xc*tanh(xc*r_hat./psi_r);
    end
    
    x_hat = SIM.beta.*x_hat_ + (1-SIM.beta)*x_hat;
    psi_x = SIM.beta*(CH.Es/2 - x_hat_.^2) + (1-SIM.beta)*psi_x;
end
y = G.xc*(2*(r_hat(1:SIM.M,:)>0)-1 + 1i*(2*(r_hat(SIM.M+1:end,:)>0)-1));
