%Find_k_step
k_step = min(k_step, length(Beliefs_sort)-1-k_init); %For Algo 2
Dum = 0;

        for k = k_init:k_step:length(Beliefs_sort)-1
            
            sum_k = sum(n_adj(k:end));
            disp = n_adj(k:end)*Beliefs_sort(k:end) - sum_k*Beliefs_sort(k);
            disp_1 = disp + sum_k*(Beliefs_sort(k)-Beliefs_sort(k+1));
            
            if disp > a*sigma^2*Zbar && disp_1 <= a*sigma^2*Zbar 
                Dum = 1; kstar = k;
                break 
            
            elseif k > k_init && disp <= a*sigma^2*Zbar
                    
                for j = k-k_step+2:2:k
                    
                    sum_j = sum(n_adj(j:end));
                    disp_j = n_adj(j:end)*Beliefs_sort(j:end) - sum_j*Beliefs_sort(j);
                    disp_j1 = disp_j + sum_j*(Beliefs_sort(j)-Beliefs_sort(j+1));
                    
                    if disp_j > a*sigma^2*Zbar && disp_j1 <= a*sigma^2*Zbar 
                        Dum = 1; kstar = j; 
                        break
                    elseif disp_j <= a*sigma^2*Zbar 
                        Dum = 1; kstar = j-1;
                        break
                    end 
                    
                end

            else
                    %do nothing
            end
            
            if Dum==1
                break
            end
            
        end

