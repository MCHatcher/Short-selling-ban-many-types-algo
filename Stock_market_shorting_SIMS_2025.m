%Stock market model with short-selling constraint and endogenous shares: simulations 
%Last updated: Apr 4, 2025. Written by Michael Hatcher (m.c.hatcher@soton.ac.uk)

clear, %clc, close all; 

%------------------
%Parameter values
%------------------
H = 1E5;  %No. of types
r = 0.1; a = 1; 
betta = 5; %3, 4.5
dbar = 0.6; sigma = 1; Zbar = 0.1;  
pf = (dbar - a*sigma^2*Zbar)/r; %Fundamental price

%----------------
%Coding choices
%----------------
Iter = 1;  %Iter = 1 turns on iterative algorithm (advisable for large H).
Fixed = 0; %Fixed  = 1: Pick fixed rather than time-varying (fitness-based) population shares. 
n_iter = 4; % no. of iterations (increase for large H)
k_step = 2; %For Algo. 2: use even number
Unconstrained = 0; %Set Unconstrained = 1 to simulate without short-selling constraints. 
T = 100;  %no. of periods

%----------------------
%Preallocate matrices
%----------------------
Bind = zeros(T,1); Bind_no = Bind; x = NaN(T,1);Check1 = NaN(T,1); 
Check11 = Check1; U = NaN(H,1); Beliefs = U; D_lag2 = U;

%--------------------------
%Generate dividend shocks 
%--------------------------
%Uncomment initially to store shocks in memory
sigma_d  = 0.005;   
pd = makedist('Normal','mu',0,'sigma',sigma_d);  %Truncated normal distribution
pd_t = truncate(pd,-dbar,dbar);
rng(1), shock = random(pd_t,T,1);   

%shock = zeros(T,1);

%-------------------------------
%Initial values and predictors 
%-------------------------------
p0 = pf + 5; x0 = p0 - pf; xlag = p0 - pf; 
%n_init = 1/H*ones(1,H); 

%Disperse beliefs 
g = zeros(H,1);
g(1:H/2) = 1.05 + 0.15*rand(H/2,1); %g = 1.2*ones(H,1); 
g(H/2+1:H) = 0;  b = zeros(H,1); C = b;
b(H/2+1:H) = linspace(-0.1,0.1,H/2); C(H/2+1:H) = 1-abs(b(H/2+1:H)); 

tic

for t=1:T 
    
    %Beliefs = NaN(H,1);
    
    if t==1
        Beliefs = b + g*x0;
        n = 1/H*ones(1,H);
    elseif t==2
        Beliefs = b + g*x(t-1);
        n = 1/H*ones(1,H);
    elseif t>=3
        Beliefs = b + g*x(t-1);
        if t==3
            Dlag2 = (b + g*x0 + a*sigma^2*Zbar - (1+r)*x(t-2))/(a*sigma^2);
        else
            Dlag2 = (b + g*x(t-3) + a*sigma^2*Zbar - (1+r)*x(t-2))/(a*sigma^2);
        end
        if Bind(t-2) == 1
            Dlag2(Dlag2<0) = 0;
        end
        U = exp(betta*( (x(t-1) + a*sigma^2*Zbar + shock(t-1) - (1+r)*x(t-2))*Dlag2 - C) );
        n = transpose(U)/sum(U);
    end
       
%------------------------------    
%Trial unconstrained solution
%------------------------------
xstar = n*Beliefs/(1+r);   

        [Beliefs_sort,I] = sort(Beliefs);  
        n_adj = n(I);

if n*Beliefs - min(Beliefs) > a*sigma^2*Zbar && Unconstrained == 0 
        
        Bind(t) = 1;
        
%Sort beliefs when there are ties (uncomment to use, not essential)
        %if length(unique(Beliefs)) ~= H
        %    run Stock_market_shorting_sort_insert
        %end
    
%--------------------------------------------       
%Obtain initial guess for no. short-sellers
%--------------------------------------------
   %Demand_star = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*xstar)/(a*sigma^2);    
   k_init0 = sum(Beliefs_sort + a*sigma^2*Zbar - (1+r)*xstar<0);
   
   %run Stock_market_shorting_iterations_insert
   if Iter == 1 
        Stock_market_shorting_iterations_insert 
   elseif Iter == 0
       k_init = k_init0;
   end
   
%-----------------------------------------   
%Find the equilibrium no.of short-sellers
%-----------------------------------------
    %run Find_k
    %run Find_k_double
    run Find_k_step
    %run Find_k_DivCon2

    Bind_no(t) = kstar;   %No. of constrained types
          
       x(t) = ( n_adj(kstar+1:end)*Beliefs_sort(kstar+1:end) - sum(n_adj(1:kstar))*a*sigma^2*Zbar  ) / ( (1+r)*sum(n_adj(kstar+1:end)) );   
       
else 
       x(t) = xstar;   %Solution when SS constraints are slack or ignored           
end

%-------------------------------------------------------
%Check market clearing (uncomment for acccuracy checks)
%-------------------------------------------------------
    %D = (Beliefs + a*sigma^2*Zbar - (1+r)*x(t))/(a*sigma^2);
    D_adj = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*x(t))/(a*sigma^2);
    if Bind(t) == 1
        %D(D<0) = 0;  
        D_adj(1:kstar) = 0; 
    end
    %Check1(t) = abs(n*D - Zbar); 
    Check11(t) = abs(n_adj*D_adj - Zbar);

end

toc

%Accuracy checks
%max(Check1)
max(Check11)

sum(Bind)
max(Bind_no), min(Bind_no)

Time = 0:T;

figure(1)
subplot(1,2,1), plot(Time,[x0 x'],'k'), hold on, title('Price deviation: x_t'), xlabel('Time, t')
subplot(1,2,2), plot(Bind_no,'k'), hold on,  title('No. of constrained types'), xlabel('Time, t'), axis([1,inf,-inf,inf])        



