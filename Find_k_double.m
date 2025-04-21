%Find_k_double

Dum = 0;

        for k = k_init:2:length(Beliefs_sort)-1
            
            sum_k = sum(n_adj(k:end));
            disp = n_adj(k:end)*Beliefs_sort(k:end) - sum_k*Beliefs_sort(k);
            disp_1 = disp + sum_k*(Beliefs_sort(k)-Beliefs_sort(k+1));
            
            if disp > a*sigma^2*Zbar && disp_1 <= a*sigma^2*Zbar 
                Dum = 1; kstar = k;
                break 
            
            elseif k > k_init && disp <= a*sigma^2*Zbar
                Dum = 1; kstar = k-1;
                break 

            else
                    %do nothing
            end
            
            if Dum==1
                break
            end
            
        end

