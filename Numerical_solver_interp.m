function [strain_num, stress_num, tau_num, s_num] = Numerical_solver(x,sp_num,de_num,dx,tau_4p, s_4p,t,e_NS,s_NS,e_PS,s_PS,strain_res)
    
    strain_num = zeros(1,length(x));
    stress_num = zeros(1,length(x));
    tau_num = zeros(1,length(x));
    s_num = zeros(1,length(x));
    for i=1:length(x)
        if i==1
            strain_num(1,i) = 0; stress_num(1,i) = 0; tau_num(1,i) = 0; s_num(1,i) = 0;
        else
            strain_num(1,i) = strain_num(1,i-1) + de_num;
            s_num(1,i) = s_num(1,i-1) + dx*(strain_num(1,i) + strain_num(1,i-1))/2;
            tau_num(1,i) = tau_slip(s_num(1,i), tau_4p, s_4p);
            stress_num(1,i) = stress_num(1,i-1) + (tau_num(1,i) + tau_num(1,i-1))/(2*t)*dx;

            %interpolation solving strain from stress
            if ismember(sp_num,[1:6,25]) % non-prestrained stress-strain
                e_temp = interp1(s_NS,e_NS,stress_num(1,i));
            else % prestrained stress-strain
                e_temp = interp1(s_PS,e_PS,stress_num(1,i));
            end

            while abs(e_temp-strain_num(1,i))>strain_res
                strain_num(1,i) = (strain_num(1,i)+e_temp)/2;
                s_num(1,i) = s_num(1,i-1) + dx*(strain_num(1,i) + strain_num(1,i-1))/2;
                tau_num(1,i) = tau_slip(s_num(1,i), tau_4p, s_4p);
                stress_num(1,i) = stress_num(1,i-1) + (tau_num(1,i) + tau_num(1,i-1))/(2*t)*dx;
                %interpolation solving strain from stress
                if ismember(sp_num,[1:6,25]) % non-prestrained stress-strain
                    e_temp = interp1(s_NS,e_NS,stress_num(1,i));
                else % prestrained stress-strain
                    e_temp = interp1(s_PS,e_PS,stress_num(1,i));
                end
                
                %monitoring the process
                residual = e_temp-strain_num(1,i);
                X = ['Specimen ', num2str(sp_num), '. Step ', num2str(i), ' /' num2str(x(end)/dx+1) '. Residual = ', num2str(double(residual))]; 
                disp(X)
           
            end
        end

    end 
    
end
    
            