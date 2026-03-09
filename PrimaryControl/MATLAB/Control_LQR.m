clc; clear variables;
Design

H = [0 0 1 0 ;
     0 0 0 1];
nr = size(H,1);

Ad_int = [Ad        Bd           ;
          H*C*Ts    zeros(nr,nu)];
Bd_int = [Bd; zeros(nr,nu)];

Qx = blkdiag(2e-2,2e-2,10,10); Qr = eye(nr);
R = eye(nu);

Q = blkdiag(Qx,Qr);
[K,~,~] = dlqr(Ad,Bd,Qx,R);
[K_int,~,~] = dlqr(Ad_int,Bd_int,Q,R);

syms z

phi = inv( z*eye(nx) - Ad );
phi_r = inv( z*eye(nx) - ( Ad - Bd*K ) );
Hcl = C*phi_r*Bd;
Hof = K*phi*Bd;

fmin_vec = 100; fmax_vec = 20e3; fvec = logspace(log10(fmin_vec), log10(fmax_vec), 500);
sigma_max = zeros(1,length(fvec));
sigma_min = zeros(1,length(fvec));

for k = 1:length(fvec)
    zi = exp(1j*2*pi*fvec(k)*Ts);
    Hi = double(subs(Hof, z, zi));
    s = svd(Hi);
    sigma_max(k) = max(s);
    sigma_min(k) = min(s);
end

set(0,'defaultTextInterpreter','latex')
fres = 1/(2*pi*sqrt(Lf*Cf));
figure(1);
semilogx(fvec,20*log10(sigma_max),'LineWidth',2); 
hold on;
semilogx(fvec,20*log10(sigma_min),'LineWidth',2); 
xline(fres,'LineWidth',2,'LineStyle','--')
xline(fsw,'LineWidth',2,'LineStyle',':')
set(gca,'FontSize',14)
hold off; grid on; xlabel('Frequency [Hz]'); ylabel('Gain [dB]')
legend('$\bar{\sigma}(\mathbf{K}_c\mathbf{\Phi}\mathbf{B}_d)$', ...
       '$\underline{\sigma}(\mathbf{K}_c\mathbf{\Phi}\mathbf{B}_d)$', ...
       'FontSize',14)