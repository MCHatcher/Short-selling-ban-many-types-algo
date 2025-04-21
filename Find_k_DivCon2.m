%Find_k_DivCon2

Dum = 0; Count = 0; 

%Initial values
k_l = k_init;
k_u = length(Beliefs_sort)-1;
k = round( k_l + ( k_u - k_l )/2 );  %to prevent overflow: k = round( ( k_l + k_u )/2 );
    
 while Dum == 0
  
    sum_k = sum(n_adj(k:end));
    disp = n_adj(k:end)*Beliefs_sort(k:end) - sum_k*Beliefs_sort(k);
    disp_1 = disp + sum_k*(Beliefs_sort(k)-Beliefs_sort(k+1));

    if disp > a*sigma^2*Zbar && disp_1 <= a*sigma^2*Zbar 
        Dum = 1; kstar = k;
    elseif disp > a*sigma^2*Zbar && disp_1 > a*sigma^2*Zbar
        
        k_l = k;
        k = ceil( k + ( k_u - k )/2 );   %to prevent overflow: k = ceil( ( k + k_u )/2 );  
    
    elseif disp <= a*sigma^2*Zbar && disp_1 <= a*sigma^2*Zbar
        
        k_u = k;
        k = floor( k_l + (k - k_l )/2 );  %to prevent overflow: k = floor( ( k_l + k )/2 );
            
    end 
        
 end
   

