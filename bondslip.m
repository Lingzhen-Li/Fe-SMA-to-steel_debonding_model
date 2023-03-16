function [slip,tau] = bondslip(slip_temp,tau_temp,i,exclude)

    if ismember(i,exclude)
        L = length(slip_temp(1,:));
%         slip_temp(:,(L-40):L) = []; 
        slip_temp(:,1:(L-40)) = [];
%         tau_temp(:,(L-40):L) = []; 
        tau_temp(:,1:(L-40)) = [];
    end
    
    if ~ismember(i,exclude)
        L = length(slip_temp(1,:));
        slip_temp(:,(L-20):L) = []; slip_temp(:,1:(L-100)) = [];
        tau_temp(:,(L-20):L) = []; tau_temp(:,1:(L-100)) = [];
    end
    tau1 = tau_temp(:);
    [slip, I] = sort(slip_temp(:)); %order slip and the position
    tau = tau1(I); % tau before smoothing
    
end