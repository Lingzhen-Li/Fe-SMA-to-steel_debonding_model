function tau = tau_slip(s, tt, ss) 
%tt is the shear stress array, constructing trilinear model
%ss is the slip array, constructing trilinear model

    % tau as a function of slip
    if s<=ss(2) % first branch
        tau = tt(2)/ss(2) * s;
    elseif s<=ss(3) % second branch
        tau = (tt(2) - tt(3))/(ss(2)-ss(3))*s + (tt(3)*ss(2) - tt(2)*ss(3))/(ss(2) - ss(3));
    elseif s<=ss(4) % third branch
        tau = tt(3)/(ss(3) - ss(4))*s - tt(3)*ss(4)/(ss(3)-ss(4));
    else
        tau = 0;
    end
    
end