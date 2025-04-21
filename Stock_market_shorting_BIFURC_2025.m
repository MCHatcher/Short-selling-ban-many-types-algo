%Stock market model with short-selling constraint and endogenous shares: bifurcation diagrams
%This version does not use a parfor loop
%Last updated: Apr 3, 2025. Written by Michael Hatcher (m.c.hatcher@soton.ac.uk).

clear; %clc; close all; 

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
k_step = 4; %For Algo. 2: use even number
Unconstrained = 0; %Set Unconstrained = 1 to simulate without short-selling constraints. 

%----------------------
%Specify belief types
%----------------------
%Parameter values
M = 4;  %No. of initial values 
window = 40;  %Sample to plot  
T = 510 + window;
num_betta = 37; num_init = 8; 
betta_min = 1; betta_max = 5; betta_split = 4.11;

%Disperse beliefs 
rng(1);
g = zeros(H,1);
g(1:H/2) = 1.1 + 0.2*rand(H/2,1); %g = 1.2*ones(H,1); 
g(H/2+1:H) = 0;  b = zeros(H,1); C = b;
b(H/2+1:H) = linspace(-0.1,0,H/2); C(H/2+1:H) = 1-abs(b(H/2+1:H)); 

%Bifurcation parameter
betta_stack1 = linspace(betta_min,betta_split,num_init); betta_stack2 = linspace(betta_split+0.005,betta_max,num_betta);
betta_stack = [betta_stack1 betta_stack2];
%betta_stack = linspace(betta_min,betta_max,num_betta);
num_betta = length(betta_stack);  dev = NaN(num_betta,1); dev1 = dev; dum = dev;

%-----------------
%Dividend shocks
%-----------------
shock = zeros(T,1);  %Deterministic skeleton

%------------------------
%Specify initial values
%------------------------
n_init = 1/H*ones(1,H);

%Baseline case
rng(3); init_stack = pf - rand(M,1);

%----------------------
%Preallocate matrices
%----------------------
brk = zeros(M,1); percent = zeros(length(betta_stack),1); C1 = brk; C11 = brk; C12 = brk; sd_x = brk; 
x_stack = NaN(window,M); x_plot = NaN(window*M,1); sd_plot1 = zeros(M,num_betta);
U = NaN(H,1); Bind = zeros(T,1); Bind_no = NaN(T,1); D = NaN(H,1); Check1 = D; Check11 = D;

tic

for v=1:num_betta 

%for v=[9 10 11 12]
    
    betta = betta_stack(v);
    x = NaN(T,1);
    
for m = 1:M

%Initial price    
p0 = init_stack(m); x0 = p0 - pf; xlag = p0 - pf;
Bind = zeros(T,1);

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
    run Find_k
    %run Find_k_double
    %run Find_k_step
    %run Find_k_DivCon2

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
        %D_adj(D_adj<0) = 0;
        D_adj(1:kstar) = 0; 
    end
    Check1(t) = abs(n*D - Zbar); 
    Check11(t) = abs(n_adj*D_adj - Zbar);

end

%Store value for bifurc diagram
x_stack(1:window,m) = x(end+1-window:end); 
%Check for no attractor
r1 = 1-isreal(x(end)); r2 = isnan(x(end)); r3 = isinf(x(end));

%-------------------------------
%Record sims with no attractor
%-------------------------------
if (r1+r2+r3)>0
    brk(m) = 1;  
end

C1(m) = max(Check1);
C11(m) = max(Check11);
C12(m) = max(Bind);
%sd_x(m) = std(x(end+1-50:end));

end

x_plot(:,v) = reshape(x_stack,1,[]);
%sd_plot1(:,v) = sd_x;
%sd_plot1(:,v) = sd_plot1(:,v)./sd_plot(:,v);
%Uncomment after simulating Unconstrained = 0 to get relative SD for plot

%Percentage of sims with no attractor
percent(v) = 100*sum(brk)/M;

dev(v) = max(C1);
dev1(v) = max(C11);
dum(v) = max(C12);

end

%Accuracy checks
%max(dev)
max(dev1)
%Check whether SS constraint binds in one or more sims
max(dum)

%-----------------
% Plot results
%-----------------

figure(1)
subplot(1,2,2)
hold on,
xlabel('Intensity of choice \beta'), ylabel('Price deviation \it{x}_t'), %title('Absence of short-selling constraints')
axis([betta_min,betta_max,-inf,inf]), set(gca, 'box','on')
plot(betta_stack, x_plot,'o', 'MarkerSize', 2.2, 'color',[0.5,0.5,0.5]) %[0.5,0.5,0.5], 'k'

toc

