clear
clc

% set plot format in accordance with latex
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');

% high quality of eps figure
set(gcf,'renderer','painters')

% default global font size
FT=14;
set(0,'defaultAxesFontSize',FT)

%% Load test data
load SMA_lap-shear_test.mat

%% Non-prestrained SMA
load Non_prestrained_SMA.mat
% read stress-strain up to 10% strain, for fitting
strain_NS = Non_prestrained_SMA(1:3000,1)/100;
stress_NS = Non_prestrained_SMA(1:3000,2);

%%
t = 1.5; % SMA thickness, in mm
b_s = 50; % SMA width, in mm

de = 0.01/100; %strain interval

a(1) = 1.777e+05; b(1) = 0.01303;
a(2) = -569.3; b(2) = -342.8;
e_NS = 0:de:20/100; %engineering strain of SMA
s_NS = a(1)*exp(b(1)*e_NS) + a(2)*exp(b(2)*e_NS)-a(1)-a(2);
ds_de_NS = a(1)*b(1)*exp(b(1)*e_NS) + a(2)*b(2)*exp(b(2)*e_NS); %tangent modulus
Gf_NS = t*(a(1)*exp(b(1)*e_NS).*(b(1)*e_NS-1)/b(1) + a(2)*exp(b(2)*e_NS).*(b(2)*e_NS-1)/b(2) + a(1)/b(1) + a(2)/b(2));

F_NS = s_NS*b_s*t/1000; % bond capacity, in kN

F_prediction([1:6 25]) = interp1(Gf_NS,F_NS,Gf_test([1:6 25])); ...
    %prediction of Fb using the experimentally measured Gf


%% Prestrained SMA
load Prestrained_SMA_Sizhe.mat
% read stress-strain up to 10% strain, for fitting
strain_PS = Prestrained_SMA(1:3000,1)/100;
stress_PS = Prestrained_SMA(1:3000,2);

%%

a(3) = 7.323e+04; b(3) = 0.02753;
a(4) = -660.2; b(4) = -239.1;
e_PS = 0:de:20/100; %engineering strain of SMA
s_PS = a(3)*exp(b(3)*e_PS) + a(4)*exp(b(4)*e_PS)-a(3)-a(4);
ds_de_PS = a(3)*b(3)*exp(b(3)*e_PS) + a(4)*b(4)*exp(b(4)*e_PS); %tangent modulus
Gf_PS = t*(a(3)*exp(b(3)*e_PS).*(b(3)*e_PS-1)/b(3) + a(4)*exp(b(4)*e_PS).*(b(4)*e_PS-1)/b(4) + a(3)/b(3) + a(4)/b(4));

F_PS = s_PS*b_s*t/1000; % bond capacity, in kN

F_prediction([7:24 26]) = interp1(Gf_PS,F_PS,Gf_test([7:24 26])); ...
    %prediction of Fb using the experimentally measured Gf

%%
figure(1)
plot(e_NS,s_NS,'b--')
hold on
plot(strain_NS,stress_NS,'b')
hold on
plot(e_PS,s_PS,'r--')
hold on
plot(strain_PS,stress_PS,'r')
hold on
xlabel('Engineering strain')
ylabel('Engineering stress (MPa)')
xlim([0 0.1])
ylim([0 900])
legend('Non-prestrained, model','Non-prestrained, tested',...
    'Prestrained, model','Prestrained, tested','Location','Southeast')
% title('$\sigma-\varepsilon$')
box on
grid on
set(gca,'GridLineStyle','--')
% grid minor
hold off

figure(2)
plot(Gf_NS, F_NS,'b--') % analytical prediction, non-prestrained
hold on
scatter(Gf_test([1:6 25]),F_test([1:6 25]),'bs') % test data, non-prestrained
hold on
plot(Gf_PS, F_PS,'r-') % analytical prediction, prestrained
hold on
scatter(Gf_test([7:24 26]),F_test([7:24 26]),'ro') % test data, prestrained
xlim([0 25])
ylim([0 70])
xlabel('Interfacial fracture energy (MPa$\cdot$mm)')
ylabel('Bond capacity (kN)')
legend('Analytical prediction, non-prestrained','Experiment, non-prestrained',...
    'Analytical prediction, prestrained','Experiment, prestrained','Location','Southeast')
box on
grid on
set(gca,'GridLineStyle','--')
% grid minor
hold off

fig_num = 2;


%% Processing load-displacement to/from three-stage bond-slip
[BSinfo, BSinfo2, BSname, ~,~] = Bond_slip_info;
load All_data_new_sigma_epsilon.mat


for sp_num=[20] % enter 11 or 20. Due to an uploading limitation, only two specimens are shown here.
    % sp_num = 20 with nonlinear adhesive, PS-S2-T0.5-2
    % sp_num = 11 with linear adhesive, PS-S1-T2-1
    if ismember(sp_num,[1:6 25]) %specimens 1-6
        Gff = Gf_NS; FF = F_NS;
    else %specimens 7-24
        Gff = Gf_PS; FF = F_PS;
    end
      
    
    %processing bond-slip from load-displacement
    s_exp = specimen{2,1,sp_num};
    F_exp = specimen{2,2,sp_num};
    Gf_new = []; tau_new = []; s_new = [];
    Gf_new = interp1(FF,Gff,F_exp,'linear','extrap');
    x_test = []; xx_test = []; strain_test = []; stress_test = []; shear_test = [];
    slip_test2 = []; stress_test2 = []; F_test2 = [];
    x_test = specimen{6,1,sp_num}; %coordinate corresponding to test SMA strain/stress
    xx_test = specimen{8,1,sp_num}; %coordinate corresponding to test shear stress
    strain_test = specimen{6,2,sp_num}; %experimental SMA tensile strain
    stress_test = specimen{7,2,sp_num}; %experimental SMA tensile stress
    shear_test = specimen{8,2,sp_num}; %experimental shear stress
    slip_test = specimen{9,2,sp_num}; %experimental slip along the bond line
    slip_test_temp = slip_test(:,50:end); % info close to the free end are useless
    stress_test_temp = stress_test(:,50:end); % info close to the free end are useless
    [slip_test2, I] = sort(slip_test_temp(:)); %re-order slip in a array
    stress_test2 = stress_test_temp(I); %re-order stress
    F_test2 = stress_test2*b_s*t/1000; % tensile force in SMA along the bond line, in kN
    
    %interval when processing data
        di=10;
    for i=1:length(Gf_new)-di
        tau_new(i) = (Gf_new(i+di)-Gf_new(i))/(s_exp(i+di)-s_exp(i));
        s_new(i) = (s_exp(i+di)+s_exp(i))/2;
        tau_new = smoothdata(tau_new,'loess',20); % window size of 3 and 10 are good.
        %do not use kernal average, because the initial point should go
        %through origin and Gauss does not.
    end
    

    % Bond-slip info from DIC,trilinear
    % tau_s_DIC, four stress points forming a bond-slip curve
    % s_DIC, four slip points forming a bond-slip curve
    tau_s_DIC(1,sp_num) = 0; s_DIC(1,sp_num) = 0; tau_s_DIC(4,sp_num) = 0;
    s_DIC(2,sp_num) = BSinfo(sp_num,1); %slip at the maximum shear, in mm 
    tau_s_DIC(2,sp_num) = BSinfo(sp_num,4); %maximum shear stress, in MPa 
    Gf1_DIC(sp_num) = tau_s_DIC(2,sp_num)*s_DIC(2,sp_num)/2;
    tau_s_DIC(3,sp_num) = BSinfo(sp_num,5); %shear stress at the turning point, in MPa
    if ismember(sp_num,[1 2 7 8 9 10 11 12])
        Gf2_DIC(sp_num) = Gf_test(sp_num);
        s_DIC(3,sp_num) = 2*Gf2_DIC(sp_num)/tau_s_DIC(2,sp_num);
        s_DIC(4,sp_num) = s_DIC(3,sp_num); %maximum slip, in mm
        Gf3_DIC(sp_num) = Gf2_DIC(sp_num) + tau_s_DIC(3,sp_num)*(s_DIC(4,sp_num)+s_DIC(3,sp_num))/2;
    else
        s_DIC(3,sp_num) = BSinfo(sp_num,2); %slip at the turning point, in mm 
        Gf2_DIC(sp_num) = Gf1_DIC(sp_num) - (tau_s_DIC(2,sp_num) - tau_s_DIC(3,sp_num))*(s_DIC(3,sp_num)+s_DIC(2,sp_num))/2 ...
            - (tau_s_DIC(3,sp_num)*s_DIC(2,sp_num) - tau_s_DIC(2,sp_num)*s_DIC(3,sp_num));
        Gf3_DIC(sp_num) = Gf_test(sp_num);
        s_DIC(4,sp_num) = 2*(Gf3_DIC(sp_num)-Gf2_DIC(sp_num))/tau_s_DIC(3,sp_num)+s_DIC(3,sp_num); %maximum slip, in mm
    end
    
    
    ss_DIC(1,:,sp_num) = linspace(0,s_DIC(4,sp_num),500); % create a space for slip, drawing bond-slip and F-s
    for i=1:length(ss_DIC(1,:,sp_num))
        if ss_DIC(1,i,sp_num)<=s_DIC(2,sp_num)
            Gf_DIC(1,i,sp_num) = Gf1_DIC(sp_num)*(ss_DIC(1,i,sp_num)/s_DIC(2,sp_num))^2;
        elseif ss_DIC(1,i,sp_num)<=s_DIC(3,sp_num)
            Gf_DIC(1,i,sp_num) = Gf1_DIC(sp_num) + (2*tau_s_DIC(2,sp_num) - (s_DIC(2,sp_num)-ss_DIC(1,i,sp_num))*...
                (tau_s_DIC(2,sp_num)-tau_s_DIC(3,sp_num))/(s_DIC(2,sp_num)-s_DIC(3,sp_num)))*(ss_DIC(1,i,sp_num)-s_DIC(2,sp_num))/2;
        else
            Gf_DIC(1,i,sp_num) = Gf2_DIC(sp_num) + (2*tau_s_DIC(3,sp_num) - (s_DIC(3,sp_num)-ss_DIC(1,i,sp_num))*...
                (tau_s_DIC(3,sp_num))/(s_DIC(3,sp_num)-s_DIC(4,sp_num)))*(ss_DIC(1,i,sp_num)-s_DIC(3,sp_num))/2;
        end
    end
    
    F_DIC(1,:,sp_num) = interp1(Gff,FF,Gf_DIC(1,:,sp_num));
       
    
    % Bond-slip info from load-displacement, trapezoidal   
    tau_s_LLZ(1,sp_num) = 0; s_LLZ(1,sp_num) = 0; tau_s_LLZ(4,sp_num) = 0;
    s_LLZ(2,sp_num) = BSinfo2(sp_num,1); %slip at the maximum shear, in mm 
    tau_s_LLZ(2,sp_num) = BSinfo2(sp_num,4); %maximum shear stress, in MPa 
    Gf1_LLZ(sp_num) = tau_s_LLZ(2,sp_num)*s_LLZ(2,sp_num)/2;
    tau_s_LLZ(3,sp_num) = BSinfo2(sp_num,5); %shear stress at the turning point, in MPa
    if ismember(sp_num,[1 2 7 8 9 10 11 12])
        Gf2_LLZ(sp_num) = Gf_test(sp_num);
        s_LLZ(3,sp_num) = 2*Gf2_LLZ(sp_num)/tau_s_LLZ(2,sp_num);
        s_LLZ(4,sp_num) = s_LLZ(3,sp_num)+1; %maximum slip, in mm. +1 or +n does not influence the result.
        Gf3_LLZ(sp_num) = Gf2_LLZ(sp_num) + tau_s_LLZ(3,sp_num)*(s_LLZ(4,sp_num)+s_LLZ(3,sp_num))/2;
    else
        s_LLZ(3,sp_num) = BSinfo2(sp_num,2); %slip at the turning point, in mm 
        Gf2_LLZ(sp_num) = Gf1_LLZ(sp_num) - (tau_s_LLZ(2,sp_num) - tau_s_LLZ(3,sp_num))*(s_LLZ(3,sp_num)+s_LLZ(2,sp_num))/2 ...
            - (tau_s_LLZ(3,sp_num)*s_LLZ(2,sp_num) - tau_s_LLZ(2,sp_num)*s_LLZ(3,sp_num));
        Gf3_LLZ(sp_num) = interp1(FF,Gff,F_test(sp_num));
        s_LLZ(4,sp_num) = 2*(Gf3_LLZ(sp_num)-Gf2_LLZ(sp_num))/tau_s_LLZ(3,sp_num)+s_LLZ(3,sp_num); %maximum slip, in mm
%         s_LLZ(4,sp_num) = BSinfo2(sp_num,3); %maximum slip, in mm
    end

    % Gf,i as a function of slip
    ss_LLZ(1,:,sp_num) = linspace(0,s_LLZ(4,sp_num),500); % create a space for slip, drawing bond-slip and F-s
    for i=1:length(ss_LLZ(1,:,sp_num))
        if ss_LLZ(1,i,sp_num)<=s_LLZ(2,sp_num)
            Gf_LLZ(1,i,sp_num) = Gf1_LLZ(sp_num)*(ss_LLZ(1,i,sp_num)/s_LLZ(2,sp_num))^2;
        elseif ss_LLZ(1,i,sp_num)<=s_LLZ(3,sp_num)
            Gf_LLZ(1,i,sp_num) = Gf1_LLZ(sp_num) + (2*tau_s_LLZ(2,sp_num) - (s_LLZ(2,sp_num)-ss_LLZ(1,i,sp_num))*...
                (tau_s_LLZ(2,sp_num)-tau_s_LLZ(3,sp_num))/(s_LLZ(2,sp_num)-s_LLZ(3,sp_num)))*(ss_LLZ(1,i,sp_num)-s_LLZ(2,sp_num))/2;
%             Gf_LLZ(1,i,sp_num) = Gf1_LLZ(sp_num) + tau_s_LLZ(2,sp_num)*(ss_LLZ(1,i,sp_num)-s_LLZ(2,sp_num))+...
%                 (ss_LLZ(1,i,sp_num)-s_LLZ(2,sp_num))^2/(s_LLZ(3,sp_num)-s_LLZ(2,sp_num))*(tau_s_LLZ(3,sp_num)-tau_s_LLZ(2,sp_num))/2;
        else
            Gf_LLZ(1,i,sp_num) = Gf2_LLZ(sp_num) + (2*tau_s_LLZ(3,sp_num) - (s_LLZ(3,sp_num)-ss_LLZ(1,i,sp_num))*...
                (tau_s_LLZ(3,sp_num))/(s_LLZ(3,sp_num)-s_LLZ(4,sp_num)))*(ss_LLZ(1,i,sp_num)-s_LLZ(3,sp_num))/2;
        end
    end
    
    F_LLZ(1,:,sp_num) = interp1(Gff,FF,Gf_LLZ(1,:,sp_num));
    
    
% Numerical solution of the lap-shear behavior
    if ismember(sp_num,[1 2 7:12]) % linear adhesive
        dx = 3; % coordinate incremental, in mm. For precise result use 3.
        x(1,:,sp_num) = 0:dx:450; % create linear space for coordinate, L=450 for linear adhesive, L=750 for nonlinear adhesive
    else % nonlinear adhesive
        dx = 5; % coordinate incremental, in mm. For precise result use 5.
        x(1,:,sp_num) = 0:dx:750; % create linear space for coordinate, L=450 for linear adhesive, L=750 for nonlinear adhesive
    end
    de_num = 1e-3; % arbitrary guess of strain incremental, in m/m. 1e-4 ~ 1e-2 are fine, did try other values
    strain_res = 1e-6;
    
    
    %Interpolation solve the lap-shear problem
    %with trapezoidal bond-slip
    [strain_num_LLZ(1,:,sp_num), stress_num_LLZ(1,:,sp_num), tau_num_LLZ(1,:,sp_num), s_num_LLZ(1,:,sp_num)] = ...
        Numerical_solver_interp(x(1,:,sp_num),sp_num,de_num,dx,tau_s_LLZ(:,sp_num), s_LLZ(:,sp_num),t,e_NS,s_NS,e_PS,s_PS,strain_res);      

    %with trilinear bond-slip
    [strain_num_DIC(1,:,sp_num), stress_num_DIC(1,:,sp_num), tau_num_DIC(1,:,sp_num), s_num_DIC(1,:,sp_num)] = ...
        Numerical_solver_interp(x(1,:,sp_num),sp_num,de_num,dx,tau_s_DIC(:,sp_num), s_DIC(:,sp_num),t,e_NS,s_NS,e_PS,s_PS,strain_res);  

    
    
    
    figure(fig_num+6*sp_num-5)
    plot(s_new,tau_new,'r') %bond-slip from load displacement
    hold on 
    exclude = [1];
    [slip_DIC, tau_DIC] = bondslip(specimen{4,1,sp_num},specimen{4,2,sp_num},sp_num,exclude);
    tau_DIC = smoothdata(tau_DIC,'loess',50); %50 is a good window size
    plot(slip_DIC,tau_DIC,'k') %bond-slip from DIC processing
    hold on
    if max(ss_DIC(1,:,sp_num))<1
        xlim([0 0.5])
    else
        xlim([0 ceil(max(ss_DIC(1,:,sp_num)))])
    end
    ylim([0 40])
    xlabel('Slip (mm)')
    ylabel('Shear stress (MPa)')
    legend('From $F-\Delta$ curve','From traditional process')
%     title('Specimen: '+ specimen{1,1,sp_num}')
    box on
    hold off
    
    grid on
    set(gca,'GridLineStyle','--')
    hold off
    
    
    figure(fig_num+6*sp_num-4)  
    subplot(2,1,1)
    plot(s_new,tau_new,'r') %bond-slip from load displacement
    hold on
    plot(s_LLZ(:,sp_num),tau_s_LLZ(:,sp_num),'b--') %simplification of bond-slip from load displacement
    hold on
    exclude = [1];
    [slip_DIC, tau_DIC] = bondslip(specimen{4,1,sp_num},specimen{4,2,sp_num},sp_num,exclude);
    tau_DIC = smoothdata(tau_DIC,'loess',50); %50 is a good window size
    if max(ss_DIC(1,:,sp_num))<1
        xlim([0 0.5])
    else
        xlim([0 ceil(max(ss_DIC(1,:,sp_num)))])
    end
    ylim([0 30])
    xlabel('Slip (mm)')
    ylabel('Shear stress (MPa)')
%     legend('From $F-\Delta$ curve','Simplified triangle')
    legend('From $F-\Delta$ curve','Simplified trapezoid')
    box on
    grid on
    set(gca,'GridLineStyle','--')
    hold off
    
    
    subplot(2,1,2)
    [s_exp, I] = sort(s_exp);
    F_exp = smoothdata(F_exp(I),'loess',50);
    plot(s_exp,F_exp,'r') %experimental load-displacement
    hold on
    plot(ss_LLZ(1,:,sp_num), F_LLZ(1,:,sp_num),'b--') %load-displacement from my model, trapezoidal bond-slip
    hold on
    xlabel('Displacement at the loaded end (mm)')
    ylabel('Force (kN)')
    if max(ss_DIC(1,:,sp_num))<1
        xlim([0 0.5])
    else
        xlim([0 ceil(max(ss_DIC(1,:,sp_num)))])
    end
    ylim([0 70])
    legend('Experiment','Modelling','Location','Southeast')
    box on
    grid on
    set(gca,'GridLineStyle','--')
    hold off   

end

%% Full-range behavior, trapezoidal, sp20
p = [23, 145, 300, 380];

if sp_num==20
    for ii=1:length(x(:,:,sp_num))
        [~,index_x_left] = min(abs(x(:,:,sp_num)-(x(:,ii,sp_num)-300)));
        if(x(:,ii,sp_num)<=300)
            F_num_LLZ(1,ii,sp_num) = stress_num_LLZ(1,ii,sp_num)*b_s*t/1000;
            strain_num_LLZ_corr(1,:,sp_num) = strain_num_LLZ(1,:,sp_num);
            s_num_LLZ_snapback(1,ii,sp_num) = s_num_LLZ(1,ii,sp_num);
        else
            F_num_LLZ(1,ii,sp_num) = trapz(x(:,index_x_left:ii,sp_num),tau_num_LLZ(1,index_x_left:ii,sp_num))*b_s/1e3;
            
            strain_num_LLZ_corr(1,index_x_left:ii,sp_num) = strain_num_LLZ(1,index_x_left:ii,sp_num)-stress_num_LLZ(1,index_x_left,sp_num)/180000;
            
            Gf_temp = interp1(FF,Gff,stress_num_LLZ(1,index_x_left,sp_num)*b_s*t/1e3);
            s_temp = interp1(Gf_LLZ(1,:,sp_num),ss_LLZ(1,:,sp_num),Gf_temp);
            s_num_LLZ_snapback(1,ii,sp_num) = s_temp + trapz(x(:,index_x_left:ii,sp_num),strain_num_LLZ_corr(1,index_x_left:ii,sp_num));
        end
    end

for i=1:length(p)
    if ismember(i,[1 2])
        [~,frame_LLZ(i)] = min(abs(F_exp(p(i)) - F_num_LLZ(1,1:length(F_num_LLZ)/2,sp_num)));
    else
        [~,frame_LLZ(i)] = min(abs(s_exp(p(i)) - s_num_LLZ_snapback(1,:,sp_num)));
    end
end
 
frame_LLZ(4) = frame_LLZ(4) - 2;

    frame_LLZ(length(p)+1)=[108]; % no later than 120 with step length of 5

        
        figure(101) % load displacement, total
        plot(s_exp, F_exp(1:length(s_exp)), 'k'); %test results
        hold on
        plot(s_num_LLZ_snapback(1,:,sp_num),F_num_LLZ(1,:,sp_num),'k--')
        hold on
        plot(s_exp(p), F_exp(p),'ko');
        hold on
        scatter(s_num_LLZ_snapback(1,frame_LLZ,sp_num),F_num_LLZ(1,frame_LLZ,sp_num),'ks')
        hold on
        xlabel('Displacement of the loaded end (mm)');
        ylabel('Load (kN)');
        xlim([0 30]); ylim([0 70]);
%         title('Specimen: '+ specimen{1,1,sp_num}')
        text(0.5,13,'(i)','FontSize',FT)
        text(1,57,'(ii)','FontSize',FT)
        text(5.3,65,'(iii)','FontSize',FT)
        text(11.5,65,'(iv)','FontSize',FT)
        text(23,23,'(v)','FontSize',FT)
        box on
        grid on
        set(gca,'GridLineStyle','--')
        hold off
        
        
        figure(102) %shear stress profile
        for i = 1:length(frame_LLZ)
            ii = frame_LLZ(i);
            [~,index_x_left] = min(abs(x(:,:,sp_num)-(x(:,ii,sp_num)-300)));
            if(x(:,ii,sp_num)<=300)
                F_num_LLZ(1,ii,sp_num) = stress_num_LLZ(1,ii,sp_num)*b_s*t/1000;
                s_num_LLZ_snapback(1,ii,sp_num) = s_num_LLZ(1,ii,sp_num);
            else
                F_num_LLZ(1,ii,sp_num) = trapz(x(:,index_x_left:ii,sp_num),tau_num_LLZ(1,index_x_left:ii,sp_num))*b_s/1e3;

                strain_num_LLZ_corr(1,index_x_left:ii,sp_num) = strain_num_LLZ(1,index_x_left:ii,sp_num)-stress_num_LLZ(1,index_x_left,sp_num)/180000;

                Gf_temp = interp1(FF,Gff,stress_num_LLZ(1,index_x_left,sp_num)*b_s*t/1e3);
                s_temp = interp1(Gf_LLZ(1,:,sp_num),ss_LLZ(1,:,sp_num),Gf_temp);
                s_num_LLZ_snapback(1,ii,sp_num) = s_temp + trapz(x(:,index_x_left:ii,sp_num),strain_num_LLZ_corr(1,index_x_left:ii,sp_num));
            end
            
            if i<=length(p)
                slip_test_temp = []; 
                slip_test_temp = slip_test(p(i),:);
%                 slip_test_back(i,:) = slip_test_temp;
                slip_test_back(i,:) = smoothdata(slip_test_temp,'gauss',20);
%                 slip_test_back(i,:) = slip_test_back(i,:)-min(slip_test_back(i,:));
                slip_test_back(i,:) = slip_test_back(i,:)-slip_test_back(i,1);
                shear_exp_back(i,:) = interp1([s_LLZ(:,sp_num); 10],[tau_s_LLZ(:,sp_num); 0],slip_test_back(i,:),'linear','extrap');
                % integral tensile stress
                shear_exp_back_temp = 1/2*([0 shear_exp_back(i,:)]+[shear_exp_back(i,:) 0]);
                shear_exp_back_temp(end) = [];
                sigma_exp_back(i,:) = (x_test(2)-x_test(1))/t* cumsum(shear_exp_back(i,:));
                % read tensile strain
                if sp_num<=6 % non-prestrained stress-strain
                    epsilon_exp_back(i,:) = interp1(s_NS,e_NS,sigma_exp_back(i,:),'linear','extrap');
                else % prestrained stress-strain
                    epsilon_exp_back(i,:) = interp1(s_PS,e_PS,sigma_exp_back(i,:),'linear','extrap');
                end
                delta_s_exp_back(i,:) = epsilon_exp_back(i,:)*(x_test(2)-x_test(1));
                s_exp_new(i,:) = cumsum(delta_s_exp_back(i,:));

                plot(x_test, shear_exp_back(i,:),'k') % test backward
                hold on
            end
            plot(x(1,:,sp_num)+300-x(:,ii,sp_num),tau_num_LLZ(1,:,sp_num),'k--') %with trapezoidal bond-slip
            hold on
        end
        xlim([0 300])
        ylim([0 30])
        xlabel('Distance to the free end (mm)')
        ylabel('Shear stress (MPa)')
        text(275,7,'(i)','FontSize',FT)
        text(270,22,'(ii)','FontSize',FT)
        text(210,22,'(iii)','FontSize',FT)
        text(135,22,'(iv)','FontSize',FT)
        text(10,22,'(v)','FontSize',FT)
        box on
        grid on
        set(gca,'GridLineStyle','--')
        hold off
        
        
        figure(103) % tensile stress profile
        for i = 1:length(frame_LLZ)
            if i<=length(p)
                plot(x_test, sigma_exp_back(i,:),'b') % test backward
                hold on
            end
            ii = frame_LLZ(i);
            [~,index_x_left] = min(abs(x(:,:,sp_num)-(x(:,ii,sp_num)-300)));
            plot(x(1,:,sp_num)+300-x(:,ii,sp_num),stress_num_LLZ(1,:,sp_num) - stress_num_LLZ(1,index_x_left,sp_num),'b--') %residual tensile stress
            hold on
        end
        xlim([0 300])
        ylim([0 1e3])
%         ylim([0 max(stress_num_LLZ(1,:,sp_num))])
        xlabel('Distance to the free end (mm)')
        ylabel('Tensile stress (MPa)')
        text(270,75,'(i)','FontSize',FT)
        text(280,450,'(ii)','FontSize',FT)
        text(219,450,'(iii)','FontSize',FT)
        text(115,450,'(iv)','FontSize',FT)
        text(20,310,'(v)','FontSize',FT)
        box on
        grid on
        set(gca,'GridLineStyle','--')
        hold off
        

        figure(104) %tensile strain profile
        for i = 1:length(frame_LLZ)
            if i<=length(p)
                plot(x_test, epsilon_exp_back(i,:),'r') % test backward
                hold on
            end
            ii = frame_LLZ(i);
            [~,index_x_left] = min(abs(x(:,:,sp_num)-(x(:,ii,sp_num)-300)));
            if(x(:,ii,sp_num)<=300)
                F_num_LLZ(1,ii,sp_num) = stress_num_LLZ(1,ii,sp_num)*b_s*t/1000;
                strain_num_LLZ_corr(1,:,sp_num) = strain_num_LLZ(1,:,sp_num);
                s_num_LLZ_snapback(1,ii,sp_num) = s_num_LLZ(1,ii,sp_num);
            else
                F_num_LLZ(1,ii,sp_num) = trapz(x(:,index_x_left:ii,sp_num),tau_num_LLZ(1,index_x_left:ii,sp_num))*b_s/1e3;

                strain_num_LLZ_corr(1,index_x_left:ii,sp_num) = strain_num_LLZ(1,index_x_left:ii,sp_num)-stress_num_LLZ(1,index_x_left,sp_num)/180000;

                Gf_temp = interp1(FF,Gff,stress_num_LLZ(1,index_x_left,sp_num)*b_s*t/1e3);
                s_temp = interp1(Gf_LLZ(1,:,sp_num),ss_LLZ(1,:,sp_num),Gf_temp);
                s_num_LLZ_snapback(1,ii,sp_num) = s_temp + trapz(x(:,index_x_left:ii,sp_num),strain_num_LLZ_corr(1,index_x_left:ii,sp_num));
            end
            plot(x(1,:,sp_num)+300-x(:,ii,sp_num),strain_num_LLZ_corr(1,:,sp_num),'r--') %residual strain
            hold on
        end
        
        xlim([0 300])
        ylim([0 0.1])
%         ylim([0 max(strain_num_LLZ(1,:,sp_num))])
        xlabel('Distance to the free end (mm)')
        ylabel('Tensile strain ($\varepsilon$)')
        text(285,0.004,'(i)','FontSize',FT)
        text(270,0.025,'(ii)','FontSize',FT)
        text(213,0.025,'(iii)','FontSize',FT)
        text(163,0.045,'(iv)','FontSize',FT)
        text(22,0.055,'(v)','FontSize',FT)
        box on
        grid on
        set(gca,'GridLineStyle','--')
        hold off
        

end

