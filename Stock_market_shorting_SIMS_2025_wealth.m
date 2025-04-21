%Stock market model with short-selling constraint and endogenous shares: simulations
%Simulation of welath and the wealth distribution across types
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
Bind = zeros(T,1); Bind_no = Bind; bound = Bind; x = NaN(T,1);Check1 = NaN(T,1); Wealth = Check1; Gini = Check1;   
Check11 = Check1; Zero_wealth = Check1; Ratio = Check1; U = NaN(H,1); Beliefs = U; D_lag2 = U; count = NaN; 
Demand_vec = NaN(H,T); Wealth_vec = Demand_vec;

%--------------------------
%Generate dividend shocks 
%--------------------------
%Uncomment initially to store shocks in memory
%sigma_d  = 0.005;   
%pd = makedist('Normal','mu',0,'sigma',sigma_d);  %Truncated normal distribution
%pd_t = truncate(pd,-dbar,dbar);
%rng(1), shock = random(pd_t,T,1);   

shock = zeros(T,1);

%-------------------------------
%Initial values and predictors 
%-------------------------------
p0 = pf + 5; x0 = p0 - pf; xlag = p0 - pf; 
%n_init = 1/H*ones(1,H); 
Wealth_init = 50*ones(H,1);

%Disperse beliefs 
g = zeros(H,1);
g(1:H/2) = 1.05 + 0.15*rand(H/2,1); %g = 1.2*ones(H,1); 
g(H/2+1:H) = 0;  b = zeros(H,1); C = b;
b(H/2+1:H) = linspace(-0.1,0.1,H/2); C(H/2+1:H) = 1-abs(b(H/2+1:H)); 

tic

for t=1:T 
    
    %Beliefs = NaN(H,1);
    %count = NaN;
    
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
    %run Find_k_DivCon2_fast

    Bind_no(t) = kstar;   %No. of constrained types
          
       x(t) = ( n_adj(kstar+1:end)*Beliefs_sort(kstar+1:end) - sum(n_adj(1:kstar))*a*sigma^2*Zbar  ) / ( (1+r)*sum(n_adj(kstar+1:end)) );   
       
else 
       x(t) = xstar;   %Solution when SS constraints are slack or ignored           
end

%-------------------------------------------------------
%Check market clearing (uncomment for acccuracy checks)
%-------------------------------------------------------
    D = (Beliefs + a*sigma^2*Zbar - (1+r)*x(t))/(a*sigma^2);
    D_adj = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*x(t))/(a*sigma^2);
    if Bind(t) == 1
        D(D<0) = 0;  
        D_adj(1:kstar) = 0; 
    end
    Check1(t) = abs(n*D - Zbar); 
    Check11(t) = abs(n_adj*D_adj - Zbar);

    Demand_vec(:,t) = D;

%-----------------------
%Track time and wealth
%-----------------------

    if t==1
        Wealth = Wealth_init;
    else
        Wealth = (1+r)*(Wealth - (x(t-1)+pf)*Demand_vec(:,t-1) ) + ( (x(t)+pf) + dbar + shock(t) )*Demand_vec(:,t-1); 
    end

    if min(Wealth) < 0
        bound(t) = 1;
    end

    Wealth(Wealth<0) = realmin;
    count(Wealth==realmin) = 1;
    Wealth_vec(:,t) = Wealth/max(Wealth);
    Wealth_norm = sort(Wealth_vec(:,t));

    %Fast approach to calculate Gini (sorted data: https://en.wikipedia.org/wiki/Gini_coefficient)
    Gini(t) =  (2/H)*(1:H)*Wealth_norm/sum(Wealth_norm) - (1+1/H);
    
    Zero_wealth(t) = sum(count); 
    Rel_wealth = Wealth/max(Wealth);

    P = prctile( Rel_wealth, [10,90] );
    Ratio(t) = P(2)/P(1);

end

x = [x0; x]; Time = 0:T; 
max(Check1)
max(Check11)
sum(bound)

%---------------
%Plot figures
%---------------
%figure(1)
%subplot(1,3,1), plot(Time,x,'--k'), hold on, title('Price deviation: x_t'), xlabel('Time, t')
%subplot(1,3,2), plot(Gini,'--k'), hold on,  title('Gini coefficient'), xlabel('Time, t'), axis([1,inf,-inf,inf]),
%subplot(1,3,3), plot(Ratio,'--k'), hold on,  title('90:10 Ratio'), xlabel('Time, t'), axis([1,inf,-inf,inf])

figure(2)
subplot(2,2,2), hold on, histogram(Wealth_vec(:,5),80,'FaceColor',[0.5,0.5,0.5])  %[0.5,0.5,0.5]

%figure(2)
%hold on, subplot(2,2,4), plot(Time(2:end),Gini,'--k'), hold on, 
%axis([-inf,inf,-inf,inf]), title('S2: Hetero. fundamentalists, b \in [-0.2,0.2], \beta = 4.5'), ylabel('Gini coefficient'), xlabel('Time')

%figure(3)
%subplot(2,3,4), hold on, histogram(Wealth_vec(:,3),'FaceColor',[0.5,0.5,0.5]), title('Wealth distribution at t=3')      
%[0.5,0.5,0.5]

toc

%Accuracy checks
%max(Check1)
max(Check11)

sum(Bind)
max(Bind_no), min(Bind_no)

     



