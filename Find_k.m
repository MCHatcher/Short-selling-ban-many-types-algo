%Find_k
%Original algorithm to find the number of short-selling constrained, k

        for k = k_init:length(Beliefs_sort)-1
            
            sum_k = sum(n_adj(k:end));
            disp = n_adj(k:end)*Beliefs_sort(k:end) - sum_k*Beliefs_sort(k);
            disp_1 = disp + sum_k*(Beliefs_sort(k)-Beliefs_sort(k+1));
            
            %if n_adj(k:end)*Beliefs_sort(k:end) - sum(n_adj(k:end))*Beliefs_sort(k) > a*sigma^2*Zbar && n_adj(k+1:end)*Beliefs_sort(k+1:end) - sum(n_adj(k+1:end))*Beliefs_sort(k+1) <= a*sigma^2*Zbar 
            
            if disp > a*sigma^2*Zbar && disp_1 <= a*sigma^2*Zbar 
                kstar = k; Dum = 1; break
            end
        end
        