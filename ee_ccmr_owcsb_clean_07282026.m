clear all; 

% ========================================================================
% Eigenfunction expansion code for a non-linear OWC with slotted rear wall
% Features: 
%     a nonlinear OWC - iterates until beta convergence
%     a linear / parameterized slotted rear wall, which dissipated energy

% Author: Clint C. M. Reyes
% Date: July 28, 2026
% ========================================================================
%% input constants 

g = 9.81; % gravity (acc)
rho = 1000; % density of water 
rho0 = 1; % density of air 

h = 8; %2.4; %water depth (m)
N = 100; % no of evansecent modes 
b = 1.0*h; % 10 slotted barrier width 

d1 = 2.5; % 3 -he = barrier 1 depth ( at x=0 )
d2 = 2.5; % barrier 2 depth ( at x=b )

by = 1;

AIval = 0.9; % real % for single value AI 

% Slotted wall parameters 
% Suh et. al., 2006 

r0 = 1;
A2 = 0.143; % centerline to centerline dist of perforation 
b2 = 0.0143; % thickness of sb 
a2 = r0*A2; 
c2 = b2/2*(1/r0-1) + 2*A2/pi*(1 - log(4*r0) + 1/3*r0^2+ 281/180*r0^4); 
l = 2*c2;
cc = 0.6+0.4*r0^2; % for cyl piles and long waves 
cf0 = (1/r0/cc-1)^2; % cf = alpha in suh  

% PTO parameters 
ccpto = 0.61; 
r0pto = 0.0125; % for pto  - might need to ch this 
cf0pto = (1/r0pto/ccpto-1)^2; % cf = alpha in suh  

% ====== calculate AI array based on constant H/L ======
khar = 0.1:0.05:5;
Tar=2*pi*(h./(g.*khar.*tanh(khar))).^0.5; 
% ======================================================
AIar = ones(1,length(Tar)) *AIval;  
Har = AIar*2;  % wave height array 
% ======================================================
omegar = 2*pi./Tar;  % omega array
kar = khar/h; % k array
L = 2*pi./kar; % L array

nitermax = 30; % beta (PTO) max iteration 

At = zeros(1,length(Tar)); 
Ar = zeros(1,length(Tar)); 
Ct = zeros(1,length(Tar)); 
Cr = zeros(1,length(Tar)); 
Ctr = zeros(1,length(Tar)); 
Crr = zeros(1,length(Tar)); 
Atr = zeros(1,length(Tar)); 
Arr = zeros(1,length(Tar)); 
Tdm = zeros(1,length(Tar)); 
Trm = zeros(1,length(Tar)); 

cb =zeros(1,length(Tar)); 
mub=zeros(1,length(Tar)); 
Qeb=zeros(1,length(Tar)); 
phib=zeros(1,length(Tar)); 
nu=zeros(1,length(Tar)); 
p=zeros(1,length(Tar)); 
phi0=zeros(1,length(Tar)); 

beta = zeros(1,length(Tar)); % dissipation coeff f(u, cf)
dbeta = zeros(1,length(Tar)); 
nibeta = zeros(1,length(Tar)); 
nibeta_pto = zeros(1,length(Tar)); 
lambda_pto = zeros(1,length(Tar)); 
mupto=zeros(1,length(Tar));  % mupto = -omeg*V0*rho*sqrt(g/h)/rho/ca^2
Powc = zeros(1,length(Tar)); 

blamb = zeros(length(Tar),1); % x axis for plotting b/L

zar = linspace(d2,h,50); % for z integration of velocity 
u3_arz = zeros(size(zar)); 
u3_arzr= zeros(size(zar)); 


betas = zeros(1, length(Tar)); % beta - Suh; not related to our beta 
Gs = zeros(1, length(Tar)); % G - Suh 
Fh = zeros(1, length(Tar)); % horizontal force 
Fh1 = zeros(1, length(Tar)); % horizontal force x=0 
Fh2 = zeros(1, length(Tar)); % horizontal force x=w
Fh_norm = zeros(1,length(Tar)); % normalized hor force 

etah = zeros(1,length(Tar)); % eta hat - chamber mean water surface for no PTO case 
                             % diffraction problem only
etahr = zeros(1,length(Tar));

Qsar = zeros(1,length(Tar)); 
crar = zeros(1,length(Tar)); 
murar = zeros(1,length(Tar)); 

% ======= Slotted wall beta calculation ======================
% solve for initial beta using suh 2006 method 
if r0<1
for ii=1:length(Tar)
    P = l*kar(ii); 

    A3 = 8*cf0/9/pi*Har(ii)*omegar(ii)*(5+cosh(2*khar(ii)))/(2*khar(ii)+sinh(2*khar(ii))); 

    a4 = (kar(ii)/omegar(ii))^2; 
    a3 = 4*kar(ii)/omegar(ii); 
    a2p = 4+P^2; 
    a1 = 0;
    a0 = -A3^2; 
    coefvct = [a4 a3 a2p a1 a0];     % Coefficient Vector
    x = roots(coefvct); 
    betas(ii) = x(x>0); 
    test(ii)= sum(x>0);

    Gs(ii) = 1/(betas(ii)/omegar(ii) - 1j*l); 
    beta(ii) = kar(ii)/Gs(ii); % convert to our beta 
end
end

% =========== beta_pto ============
% beta_pto = beta; % assume initial guess similar to beta_sb
beta_pto = ones(size(Tar)); % initial guess for nonlinear turbine 

%% ==========================================

% eigenfunction products 
% same index 
f00a = @(a0,h) (cos(a0*h)*sin(a0*h) + a0*h)/(2*a0); 
f00b = @(a0,h,d) (cos(a0*h)*sin(a0*h) + a0*h)/(2*a0) - ... 
                 (cos(a0*(h-d))*sin(a0*(h-d)) + a0*(h-d))/(2*a0);
f00c = @(a0,h,d) (cos(a0*(h-d))*sin(a0*(h-d)) + a0*(h-d))/(2*a0);
% different index 
f01a = @(a0,a1,h)   sin((a0-a1)*h)/(2*(a0-a1)) +  sin((a0+a1)*h)/(2*(a0+a1));
f01b = @(a0,a1,h,d) sin((a0-a1)*h)/(2*(a0-a1)) +  sin((a0+a1)*h)/(2*(a0+a1)) ...
                    - ( sin((a0-a1)*(h-d))/(2*(a0-a1)) +  sin((a0+a1)*(h-d))/(2*(a0+a1)) );
f01c = @(a0,a1,h,d) sin((a0-a1)*(h-d))/(2*(a0-a1)) +  sin((a0+a1)*(h-d))/(2*(a0+a1));

% define function for cic1 - eigenfunction integrated form ( analytic and not numerical ) - for radiation problem 
f1c = @(a1,h,d) sin(a1*(h-d))/a1; % integral of eigenfunction using imag waveno 

% eigenfunctions - non integrated 
phi1 = @(k,z,h) cos(k*(z+h)); 

ii=1; % placeholder - remove later 
for kk=1:length(Tar) 

T=Tar(kk); % get wave period T from T array 
omeg = 2*pi/T;  % wave angular frequency 

% calculate wave number and kh for plotting (blamb) 
k_ar = disp_rel_ee2(T,h,g,N); % all in real form 
a0 = -1j*k_ar(1); % k0 = propagating mode wave no. in complex form
am = [a0; k_ar(2:N+1)]; % changed to incl 0 in ee_ccmr5

km = k_ar;

phi0(kk) = -1j*g/omeg/cosh(k_ar(1)*h); % wave amp to potential amp conversion factor 

blamb(kk) = k_ar(1)*h;
% blamb(kk) = k_ar(1)*b/2/pi; 

%% define eigenfunctions 
% f00 = @(z,a0,a1,h) ( cos(a0*(z+h))*sin(ao*(z+h)) + ao*(z+h) ) / (2*a0); 
% h and d are positive values 
% a0 is imag 

% n is 0-N row = eqns
% m is 1-N col = coeffs
cnma=zeros(N+1,N+1);  % m is from 0 
cnmb1=zeros(N+1,N+1);
cnmc1=zeros(N+1,N+1); 
cnmb2=zeros(N+1,N+1);
cnmc2=zeros(N+1,N+1); 
cnc1 = zeros(N+1,1); 
cnc2 = zeros(N+1,1); 

for xx=1:N+1 % n loop 
    for yy=1:N+1 % m loop 
        if xx==yy
            cnma(xx,yy) = f00a(am(xx),h); 
            cnmb1(xx,yy) = f00b(am(xx),h,d1);
            cnmc1(xx,yy) = f00c(am(xx),h,d1);
            cnmb2(xx,yy) = f00b(am(xx),h,d2); % use one d for now 
            cnmc2(xx,yy) = f00c(am(xx),h,d2);
        else
            cnma(xx,yy) = f01a(am(xx),am(yy),h); 
            cnmb1(xx,yy) = f01b(am(xx),am(yy),h,d1); 
            cnmc1(xx,yy) = f01c(am(xx),am(yy),h,d1); 
            cnmb2(xx,yy) = f01b(am(xx),am(yy),h,d2); % use one d for now 
            cnmc2(xx,yy) = f01c(am(xx),am(yy),h,d2); 
        end
    end
    cnc1(xx) = f1c(am(xx),h,d1); 
    cnc2(xx) = f1c(am(xx),h,d2); 
end

MR11 = zeros(N+1,N+1); 
MR12 = zeros(N+1,N+1); 
MR21 = zeros(N+1,N+1); 
MR22 = zeros(N+1,N+1);

MA11 = zeros(N+1,N+1); 
MA12 = zeros(N+1,N+1); 
MA21 = zeros(N+1,N+1); 
MA22 = zeros(N+1,N+1);

MB11 = zeros(N+1,N+1); 
MB12 = zeros(N+1,N+1); 
MB21 = zeros(N+1,N+1); 
MB22 = zeros(N+1,N+1);

MT11 = zeros(N+1,N+1); 
MT12 = zeros(N+1,N+1); 
MT21 = zeros(N+1,N+1); 
MT22 = zeros(N+1,N+1);

BI11 = zeros(N+1,1); 
BI12 = zeros(N+1,1); 
BI21 = zeros(N+1,1); 
BI22 = zeros(N+1,1);

BR11 = zeros(N+1,1); 
BR12 = zeros(N+1,1); 
BR21 = zeros(N+1,1); 
BR22 = zeros(N+1,1);

for xx=1:N+1
    BI11(xx) = -1j*k_ar(1)*phi0(kk)*cnma(xx,1); 
    BI12(xx) = -phi0(kk)*cnmc1(xx,1); 
    BR12(xx) = 1j/rho/omeg*cnc1(xx); % check sign later 
    BR22(xx) = -1j/rho/omeg*cnc2(xx);

    % for yy=2:N+1
    for yy=1:N+1
        % check all 22's 
        MR11(xx,yy) = am(yy)*cnma(xx,yy); 
        MR12(xx,yy) = cnmc1(xx,yy); 
        
        MA11(xx,yy) = am(yy)*cnmc1(xx,yy); 
        MA12(xx,yy) = -am(yy)*cnmb1(xx,yy) - cnmc1(xx,yy); 
        MA21(xx,yy) = -am(yy)*cnma(xx,yy)*exp(-am(yy)*b); 
        MA22(xx,yy) = cnmc2(xx,yy)*exp(-am(yy)*b); 

        MB11(xx,yy) = -am(yy)*exp(-am(yy)*b)*cnmc1(xx,yy); 
        MB12(xx,yy) = am(yy)*exp(-am(yy)*b)*cnmb1(xx,yy) - exp(-am(yy)*b)*cnmc1(xx,yy); 
        MB21(xx,yy) = am(yy)*cnma(xx,yy); 
        MB22(xx,yy) = cnmc2(xx,yy); 

        MT21(xx,yy) = am(yy)*cnmc2(xx,yy); 
        MT22(xx,yy) = -(cnmc2(xx,yy) + am(yy)*cnmb2(xx,yy) + 1j*beta(kk)/k_ar(1)*am(yy)*cnmc2(xx,yy)); 

    end
end

M = [MR11 MA11 MB11 MT11;
     MR12 MA12 MB12 MT12; 
     MR21 MA21 MB21 MT21; 
     MR22 MA22 MB22 MT22];

B = [BI11; BI12; BI21; BI22]; 

Br = -[BR11; BR12; BR21; BR22]; 

coefs = linsolve(M,B);
coefs_r = linsolve(M,Br);

At(ii,kk) = coefs(3*N+4);
Ar(ii,kk) = coefs(1); 
Atr(ii,kk) = coefs_r(3*N+4);
Arr(ii,kk) = coefs_r(1); 

Rdm = coefs(1:N+1); 
Rrm = coefs_r(1:N+1); 
Tdm = coefs(3*N+4:4*N+4); 
Trm = coefs_r(3*N+4:4*N+4); 

AL1 = coefs(N+2); 
AL1n = coefs(N+3 : 2*N+2); 
AR1 = coefs(2*N+3); 
AR1n = coefs(2*N+4 : 3*N+3); 

AL1_r = coefs_r(N+2); 
AL1n_r = coefs_r(N+3 : 2*N+2); 
AR1_r = coefs_r(2*N+3); 
AR1n_r = coefs_r(2*N+4 : 3*N+3);

% ======== phi_D and phi_R are solved ============

% calculate parameters for solving radaition problem 
Iz_ar = sin(am.*h);  % rad problem - simplified (possibly same with the scattering too! - but check! ) 
Ix_ar = exp(-am*b)-1; 
ALn_r = [AL1; AL1n]; % combine - simplified version 
ARn_r = [AR1; AR1n]; 
Qs = by * sum(Iz_ar.*Ix_ar.*(ALn_r+ARn_r)); % actually velocity, U 01/09/2025

% ==== ccmr 9/19/2024 - simplified rad problem ====
% Iz_arr = sin(am.*h);  % rad problem - simplified (possibly same with the scattering too! - but check! ) 
% Ix_arr = exp(-am*b)-1; 
ALn_rr = [AL1_r; AL1n_r]; % combine - simplified version 
ARn_rr = [AR1_r; AR1n_r]; 
Qrnop = by * sum(Iz_ar.*Ix_ar.*(ALn_rr+ARn_rr));

% obtain c and mu from radiation problem 
c = -real(Qrnop); 
mu = imag(Qrnop); 
% 
% % ---- calculate cpto ----- (new - nonlinear) 
% cpto = kar(kk)*b/rho0/omeg/beta_pto(kk)*0.01;  % compare scale in magn of old cpto and new 
% mupto = 0; % not zero in he - effect not sign 
%       cpto effect might be stronger anyway
% --------------------------
% this is for a linear PTO! rewrite after obtaining beta! 
% pto characteristics / parameters 
% mupto = omeg*V0/rho0/ca^2; 
mupto = 0; 
% calc p 
% p(kk) = Qs/(cpto+c - 1j*(mu+mupto)); % complex pressure amplitude 
% pb = cpto/2*abs(Qs)^2/((cpto+c)^2+(mu+mupto)^2); % Powc =    %% recheck! AI must be in here? 12/7/2024
Cg = omeg/2/k_ar(1)*(1+2*k_ar(1)*h/sinh(2*k_ar(1)*h)); 
Pi = 0.5*rho*g*AIar(kk)^2*Cg; 

% full velocity 
% Q = AIar(kk)*Qs + p(kk)*Qrnop; 

% nu(kk) = pb/PI; %must get actual Ai!   %% recheck! 12/7/2024

% actual Qr (*p)
% Qr = p*Qrnop; 

% cb(kk) = rho*g*c/omeg/b; 
% mub(kk) = rho*g*mu/omeg/b; 
% Qeb(kk) = abs(Qs)/omeg/b/AIar(kk); % use actual Ai!  
% phib(kk) = angle(Qs); 

% !initialize beta_pto
% beta_pto(kk) = cf0/2 * k_ar(1)*AIar(kk); 

% beta_pto iteration 
if cf0pto>0
    for ibeta=1:nitermax
    
    if ibeta==1
        % !already initialized - do nothing 
        % beta_pto(kk) = cf0/2 * k_ar(1)*AIar(kk); 
        betanew = cf0pto/2 * k_ar(1)*AIar(kk);
        % beta_pto(kk)=betanew;
    else
        betanew = (cf0pto/2)*(8/3/pi)*(k_ar(1)/omeg)*abs(U0);
    end

% calc cpto 
ce = beta_pto(kk)*omeg/kar(kk);
cpto = by * b / (rho0*ce); 
% calc Pa 
p(kk) = AIar(kk)*Qs/((cpto+c)-1j*(mupto+mu)); 
% p(kk) = 0;
Powc(kk) = cpto/2*(abs(p(kk)))^2; 
lambda_pto(kk) = Powc(kk)/Pi; 
nu(kk) = Powc(kk)/Pi/by;

% calc Q 
% Q = AIar(kk)*Qs*1 + p(kk)*Qrnop*1; % 1 = width into the page 
% calc U0 from Q
% U0 = Q/b;
U0 = AIar(kk)*Qs / ( rho0*ce * ((cpto+c)-1j*(mupto+mu)) );

dbeta(kk) = abs( betanew - beta_pto(kk) )/beta_pto(kk);

Qsar(kk) = Qs; 
crar(kk) = c; 
murar(kk) = mu; 

if dbeta(kk)<1e-3
    break
else
    beta_pto(kk) = betanew; 
    nibeta_pto(kk)=ibeta; 
end

% end % end else

end % end beta iter 

else
    beta_pto(kk)=0;
end

% ======= force calculation Fh =======
% J = sin(am(kk) 

% calc force at x = 0
Id1 = phi0(kk)*1./am(1).*sin(am(1)*h) + sum( 1./am.*sin(am*h) .* Rdm);  % check (1.0j * am(1))
% Id1 = phi0(kk)*1./(1.0j * am(1)).*sin(am(1)*h) + sum( 1./am.*sin(am*h) .* Rdm);  % check 
Id2 = sum( 1./am.*sin(am*h) .* ( ALn_r + ARn_r.*exp(-am*b) ) );
Ir1 = sum( 1./am.*sin(am*h) .* Rrm);
Ir2 = sum( 1./am.*sin(am*h) .* ( ALn_rr + ARn_rr.*exp(-am*b) ) ) + 1.0j*h/rho/omeg;

Fh1(kk) = ( 1.0j*omeg*rho* 1.0 * ( (AIar(kk))*(Id1 - Id2) + (p(kk))*(Ir1 - Ir2) ) ); % we set time phase e^-iwt = 1.0 

% calc force at x = w
Id2w = sum( 1./am.*sin(am*h) .* ( ALn_r.*exp(-am*b) + ARn_r ) );
Ir2w = sum( 1./am.*sin(am*h) .* ( ALn_rr.*exp(-am*b) + ARn_rr ) ) + 1.0j*h/rho/omeg;
Id3w = sum( 1./am.*sin(am*h) .* Tdm);
Ir3w = sum( 1./am.*sin(am*h) .* Trm);

Fh2(kk) = ( 1.0j*omeg*rho* 1.0 * ( (AIar(kk))*(Id2w - Id3w) + (p(kk))*(Ir2w - Ir3w) ) ); % we set time phase e^-iwt = 1.0 

% notes
% AI is real !
% p is not real! 

Fh(kk) = (Fh1(kk) + Fh2(kk)); 
Fh_norm(kk) = abs(Fh(kk))/ (rho*g*AIar(kk)*h); 

% calc mean water surface (wo time factor) 
etah(kk) = 1.0j / omeg / b * sum( (ALn_r + ARn_r).*(exp(-am*b) - 1).*sin(am*h) );
etahr(kk) = 1.0j / omeg / b * sum( (ALn_rr + ARn_rr).*(exp(-am*b) - 1).*sin(am*h) );

end % end kk / period loop 

% ====== end of kk looping ======= 
%%
for kk=1:length(Tar)
    for ii=1:1 
        T = Tar(kk);
        omeg = 2*pi/T; 

        Ct(ii,kk) = abs(At(ii,kk)/phi0(kk));  % scattering 
        Cr(ii,kk) = abs(Ar(ii,kk)/phi0(kk)); 

        Ctr(ii,kk) = abs( (AIar(kk)*At(ii,kk) + p(kk)*Atr(ii,kk))/ (phi0(kk)*AIar(kk)) );  % total
        Crr(ii,kk) = abs( (AIar(kk)*Ar(ii,kk) + p(kk)*Arr(ii,kk))/ (phi0(kk)*AIar(kk)) );  
    end
end
% energy balance 
checkE1 = Ct.^2 + Cr.^2 ;
checkE2 = Ctr.^2 + Crr.^2 + nu;

%%
etahn=abs(etah*AIval)/(2*AIval);
etahrn=abs(etahr.*p)/(2*AIval);
% blamb'
%% 
outdir = 'outfiles';
if ~exist(outdir)
    mkdir(outdir);
end

xcap=5;
% 
figure('position', [0 0 200 800])
subplot(6,1,1)
plot(blamb,Ctr,'bo-')
hold on 
plot(blamb,Ct,'ro-')
xlabel('kh')
ylabel('C_T')
xlim([0,xcap])
ylim([0,1])
legend('comb','scat')

subplot(6,1,2)
plot(blamb,Crr,'bo-')
hold on 
plot(blamb,Cr,'ro-')
xlabel('kh')
ylabel('C_R')
xlim([0,xcap])
ylim([0,1])
subplot(6,1,3)
hold on
plot(blamb, Fh_norm,'b-', 'linewidth',2)
xlabel('kh')
ylabel('C_D')
xlim([0,xcap])

subplot(6,1,4)
plot(blamb,abs(Qsar),'bo-')
xlabel('kh')
ylabel('$\eta$','Interpreter','Latex')
xlim([0,xcap])
% ylim([0,1])

subplot(6,1,5)
hold on 
plot(blamb,crar,'ro-') % store it and output to file 
% plot(blamb,abs(etah),'ro-')
% plot(blamb,abs(etahr),'bo-')
xlabel('kh')
ylabel("$\xi'_D$",'Interpreter','Latex')
xlim([0,xcap])
% ylim([0,1])

subplot(6,1,6)
hold on 
plot(blamb,murar,'bo-')
xlabel('kh')
ylabel("$\xi'_R$",'Interpreter','Latex')
xlim([0,xcap])
saveas(gcf, [outdir '/outallCTCR_2bar_b4_wrad_wpp2_h3mod5' num2str(d1) 'w' num2str(b) '.png'])

bL = b./L;

%% analyze peak of efficiency 
[mnu,imnu]=max(nu); % gives peak efficiency value and corr index 
Lb_peak_eta = 1/bL(imnu); % max efficiency L/b
kh_peak_eta = blamb(imnu);
%% out all to file 
fileID = fopen([outdir '/results_ee_5b4_wrad2_wpp2_h3mod6b_' 'w' num2str(b,'%04.1f') 'r0' num2str(r0,'%.2f') 'd' num2str(d1,'%.2f') 'a' num2str(AIval,'%.4f') '.txt'],'w');
fprintf(fileID, 'kh\tb/L\tCtr\tCrr\tnu\tCt\tCr\tCd\tlambda\tPowc\t |Qs| \t c \t mu\n');
% note Fh1 and Fh2 outputs only real parts! 

% note blamb is kh 
% bL is b/L 

for ii=1:length(blamb)
    fprintf(fileID,'%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\n',blamb(ii),bL(ii),Ctr(ii),Crr(ii),nu(ii),Ct(ii),Cr(ii),Fh_norm(ii),lambda_pto(ii),Powc(ii), abs(Qsar(ii)), crar(ii), murar(ii)); 
end
fclose(fileID);
