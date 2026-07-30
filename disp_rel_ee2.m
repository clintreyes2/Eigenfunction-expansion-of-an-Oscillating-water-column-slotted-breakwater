function k_ar=disp_rel_ee2(T,h,g,N)
% solves dispersion relation including evanescent (imag) wave numbers 
% N is the no of terms in evanescent mode expansion 
% T wave period 
% h water depth 
% g acc due to gravity

% k_ar are all real in output, k1-kn are actually imag

k_ar = zeros(N+1, 1); % array for wave nos. 

omeg = 2*pi/T; % angular wave freq 
% if T<0.9
k0 = 2*pi/1.56/T^2; % initial guess for waveno  
% else
%     % k0 = 
% k0 = 2*pi/sqrt(9.81*h)/T; 
% end

% newton rhapson iteration convergence conditions 
maxiter = 100000; 
eps = 1e-16; 

% solve for real k 
C = h*omeg^2/g;
f = @(x) x*tanh(x) - C; 
df = @(x) tanh(x) + x*sech(x)^2; 

x0 = k0*h; % initial guess k0*h 
% 0.9 

for i=1:maxiter 
    x1=x0-f(x0)/df(x0);
         % abs(x1-x0)
    if abs(x1-x0) < eps
        break
    else 
        x0=x1; 
    end
end
% i
k_ar(1) = x1/h; % take forward root, prev root is more standard! (shouldnt matter since tol is low) 

% solve for imag k 
% NOTE: SOLVED AS REAL K 
% k = -ik0 ( where k0 is the real k; k0 = ik ) 
% - note this is pos root. actually, k0 = +- ik  

% C = h*omeg^2/g; % no need to repeat, C is the same 
f = @(x) -x*tan(x) - C;
df = @(x) -tan(x) - x*sec(x)^2; 

% for ii=1:N
for ii=1:N
    % if ii==1
    %     x0 = -k_ar(ii)*h; 
    %     % x0 = 0.5
    % else 
    %     % x0 = imag(k_ar(ii))*h;
    % end
    % x0 = (2*ii+1)*pi/2-pi/4;  % -1 is a negative root 
    % pi/8
    % x0 = ii*pi;
    % x0 = ii*pi-pi/4;
    % last one 09-08-2024
    % x0 = (ii-0.5)*pi+pi/2; %0.5 is a factor, uses a pi period
    
    % %1 
    % if k_ar(1)>5.3%7
    % x0 = (ii-0.5+4)*pi+pi/2/2; % for short waves 
    % else
    % x0 = (ii-0.5)*pi+pi/2/2; % for long waves
    % end
    % %1 modified T 
    % % if T<0.87
    % % x0 = (ii-0.5+4)*pi+pi/2/2; % for short waves 
    % % else
    % % x0 = (ii-0.5)*pi+pi/2/2; % for long waves
    % % end
    % % x0 = (ii-0.5)*pi*2+pi*0.7; % use of a 2pi period 
    % % orig 
    % x0 = (ii-0.5)*pi+pi/4; 
    x0 = (ii-0.5)*pi*1.01; 

    % 1.0 works for 5b4
    % it depends on the slope for tan
    for i=1:maxiter
        x1=x0-f(x0)/df(x0);
    
        if abs(x1-x0) < eps
            break
        else 
            x0=x1; 
        end
    end
    % k_ar(ii+1) = 1j*x1/h; % assuming x1 is real here (but actually imag) 
    k_ar(ii+1) = x1/h; % assuming x1 is real here (but actually imag) 
    % k_ar(ii+1) = abs(x1)/h;
    % i
end
% ii