Rf = 65e-3;
Lf = 2.9e-3;
Cf = 12e-6;
Rg = 100e-3;
Lg = 3e-3;
Vg = 230*sqrt(2);
fs = 50; ws = 2*pi*fs;
fsw = 5e3; Ts = 1/(2*fsw);

A = [-Rf/Lf     ws     -1/Lf    0    ;
     -ws       -Rf/Lf   0      -1/Lf ;
      1/Cf      0       0       ws   ;
      0         1/Cf   -ws      0   ];
B = [ 1/Lf      0    ;
      0         1/Lf ;
      0         0    ;
      0         0   ];
nx = size(B,1); nu = size(B,2);
C = eye(nx);
D = zeros(nx,nu);

ss_VSI = ss(A,B,C,D); ss_VSId = c2d(ss_VSI,Ts); [Ad,Bd,~,~] = ssdata(ss_VSId);