function [BSinfo, BSinfo2, BSname, BS_pm, BS_pm_name] = Bond_slip_info()
% Each row represents a specimen
% First column: s1, in mm
% Second column: s2, in mm
% Third column: average Gf, in MPa*mm
% Fourth column: tau_max, in MPa, it is actually useless here and will be estimated from Gf 
% Fifth column: tau_t, in MPa
% Cell: sample name

% BSinfo, processed from DIC
% BSinfo2, processed from my model

%% Bond-slip from my DIC
% Non-prestrained SMA
BSinfo(1,1) = 0.03; BSinfo(1,2) = 0.1; BSinfo(1,3) = BSinfo(1,2); 
BSinfo(1,4) = 28; BSinfo(1,5) = 0; BSname{1} = 'NS-S1-T0.5-1';
BSinfo(2,1) = 0.024; BSinfo(2,2) = 0.1; BSinfo(2,3) = BSinfo(2,2); 
BSinfo(2,4) = 29; BSinfo(2,5) = 0; BSname{2} = 'NS-S1-T0.5-2';
BSinfo(3,1) = 0.065; BSinfo(3,2) = 0.12; BSinfo(3,3) = 2; 
BSinfo(3,4) = 25; BSinfo(3,5) = 10; BSname{3} = 'NS-A-T0.5-1';
BSinfo(4,1) = 0.07; BSinfo(4,2) = 0.12; BSinfo(4,3) = 2; 
BSinfo(4,4) = 26; BSinfo(4,5) = 10; BSname{4} = 'NS-A-T0.5-2';
BSinfo(5,1) = 0.06; BSinfo(5,2) = 0.11; BSinfo(5,3) = 2.5; 
BSinfo(5,4) = 31; BSinfo(5,5) = 12; BSname{5} = 'NS-S2-T0.5-1';
BSinfo(6,1) = 0.06; BSinfo(6,2) = 0.11; BSinfo(6,3) = 2.5; 
BSinfo(6,4) = 30; BSinfo(6,5) = 12; BSname{6} = 'NS-S2-T0.5-2';

% Prestrained SMA,Sika 30
BSinfo(7,1) = 0.017; BSinfo(7,2) = 0.12; BSinfo(7,3) = BSinfo(7,2); 
BSinfo(7,4) = 20; BSinfo(7,5) = 0; BSname{7} = 'PS-S1-T0.5-1';
BSinfo(8,1) = 0.018; BSinfo(8,2) = 0.12; BSinfo(8,3) = BSinfo(8,2); 
BSinfo(8,4) = 21; BSinfo(8,5) = 0; BSname{8} = 'PS-S1-T0.5-2';
BSinfo(9,1) = 0.021; BSinfo(9,2) = 0.11; BSinfo(9,3) = BSinfo(9,2); 
BSinfo(9,4) = 21; BSinfo(9,5) = 0; BSname{9} = 'PS-S1-T1-1';
BSinfo(10,1) = 0.022; BSinfo(10,2) = 0.13; BSinfo(10,3) = BSinfo(10,2); 
BSinfo(10,4) = 20; BSinfo(10,5) = 0; BSname{10} = 'PS-S1-T1-2';
BSinfo(11,1) = 0.035; BSinfo(11,2) = 0.1; BSinfo(11,3) = BSinfo(11,2); 
BSinfo(11,4) = 21; BSinfo(11,5) = 0; BSname{11} = 'PS-S1-T2-1';
BSinfo(12,1) = 0.026; BSinfo(12,2) = 0.14; BSinfo(12,3) = BSinfo(12,2); 
BSinfo(12,4) = 18; BSinfo(12,5) = 0; BSname{12} = 'PS-S1-T2-2';

% Prestrained SMA,Araldite 2015
BSinfo(13,1) = 0.15; BSinfo(13,2) = 0.23; BSinfo(13,3) = 2.5; 
BSinfo(13,4) = 32; BSinfo(13,5) = 7; BSname{13} = 'PS-A-T0.5-1';
BSinfo(14,1) = 0.15; BSinfo(14,2) = 0.23; BSinfo(14,3) = 2.5; 
BSinfo(14,4) = 35; BSinfo(14,5) = 7; BSname{14} = 'PS-A-T0.5-2';
BSinfo(15,1) = 0.21; BSinfo(15,2) = 0.30; BSinfo(15,3) = 3; 
BSinfo(15,4) = 28; BSinfo(15,5) = 9; BSname{15} = 'PS-A-T1-1';
BSinfo(16,1) = 0.21; BSinfo(16,2) = 0.30; BSinfo(16,3) = 2.5; 
BSinfo(16,4) = 28; BSinfo(16,5) = 9; BSname{16} = 'PS-A-T1-2';
BSinfo(17,1) = 0.23; BSinfo(17,2) = 0.30; BSinfo(17,3) = 3; 
BSinfo(17,4) = 26; BSinfo(17,5) = 8; BSname{17} = 'PS-A-T2-1';
BSinfo(18,1) = 0.22; BSinfo(18,2) = 0.3; BSinfo(18,3) = 2.5; 
BSinfo(18,4) = 27; BSinfo(18,5) = 8; BSname{18} = 'PS-A-T2-2';

% Prestrained SMA,Sika 1277
BSinfo(19,1) = 0.12; BSinfo(19,2) = 0.2; BSinfo(19,3) = 2.6; 
BSinfo(19,4) = 38; BSinfo(19,5) = 8; BSname{19} = 'PS-S2-T0.5-1';
BSinfo(20,1) = 0.14; BSinfo(20,2) = 0.2; BSinfo(20,3) = 2.7; 
BSinfo(20,4) = 40; BSinfo(20,5) = 9; BSname{20} = 'PS-S2-T0.5-2';
BSinfo(21,1) = 0.16; BSinfo(21,2) = 0.22; BSinfo(21,3) = 3.3; 
BSinfo(21,4) = 37; BSinfo(21,5) = 11; BSname{21} = 'PS-S2-T1-1';
BSinfo(22,1) = 0.17; BSinfo(22,2) = 0.23; BSinfo(22,3) = 3.4; 
BSinfo(22,4) = 37; BSinfo(22,5) = 11; BSname{22} = 'PS-S2-T1-2';
% BSinfo(23,1) = ; BSinfo(23,2) = ; BSinfo(23,3) = ; 
% BSinfo(23,4) = ; BSinfo(23,5) = ; BSname{23} = 'PS-S2-T2-1';
BSinfo(24,1) = 0.20; BSinfo(24,2) = 0.26; BSinfo(24,3) = 4; 
BSinfo(24,4) = 30; BSinfo(24,5) = 12; BSname{24} = 'PS-S2-T2-2';
BSinfo(23,:) = BSinfo(24,:); BSname{23} = BSname{24}; 
%'PS-S2-T2-1' did not have information, therefore, use 'PS-S2-T2-2' to represent

% Prestrained SMA,Sika 888
BSinfo(25,1) = 0.05; BSinfo(25,2) = 0.20; BSinfo(25,3) = 1.9; 
BSinfo(25,4) = 32; BSinfo(25,5) = 8.5; BSname{25} = 'NS-S3-T0.5-1';
BSinfo(26,1) = 0.14; BSinfo(26,2) = 0.35; BSinfo(26,3) = 2.3; 
BSinfo(26,4) = 33; BSinfo(26,5) = 6; BSname{26} = 'PS-S3-T0.5-1';

%% Bond-slip from my model
% Non-prestrained SMA
BSinfo2(1,1) = 0.04; BSinfo2(1,2) = 0.3; BSinfo2(1,3) = BSinfo2(1,2); 
BSinfo2(1,4) = 18; BSinfo2(1,5) = 0; BSname{1} = 'NS-S1-T0.5-1';
BSinfo2(2,1) = 0.055; BSinfo2(2,2) = 0.12; BSinfo2(2,3) = BSinfo2(2,2); 
BSinfo2(2,4) = 25; BSinfo2(2,5) = 0; BSname{2} = 'NS-S1-T0.5-2';
BSinfo2(3,1) = 0.1; BSinfo2(3,2) = 0.4; BSinfo2(3,3) = 0.8; 
BSinfo2(3,4) = 19; BSinfo2(3,5) = 19; BSname{3} = 'NS-A-T0.5-1';
BSinfo2(4,1) = 0.04; BSinfo2(4,2) = 0.5; BSinfo2(4,3) = 0.8; 
BSinfo2(4,4) = 17; BSinfo2(4,5) = 17; BSname{4} = 'NS-A-T0.5-2';
BSinfo2(5,1) = 0.04; BSinfo2(5,2) = 0.5; BSinfo2(5,3) = 1.1; 
BSinfo2(5,4) = 23; BSinfo2(5,5) = 23; BSname{5} = 'NS-S2-T0.5-1';
BSinfo2(6,1) = 0.05; BSinfo2(6,2) = 0.5; BSinfo2(6,3) = 1.1; 
BSinfo2(6,4) = 23; BSinfo2(6,5) = 23; BSname{6} = 'NS-S2-T0.5-2';

% Prestrained SMA,Sika 30
BSinfo2(7,1) = 0.03; BSinfo2(7,2) = 0.11; BSinfo2(7,3) = BSinfo2(7,2); 
BSinfo2(7,4) = 22; BSinfo2(7,5) = 0; BSname{7} = 'PS-S1-T0.5-1';
BSinfo2(8,1) = 0.04; BSinfo2(8,2) = 0.13; BSinfo2(8,3) = BSinfo2(8,2); 
BSinfo2(8,4) = 20; BSinfo2(8,5) = 0; BSname{8} = 'PS-S1-T0.5-2';
BSinfo2(9,1) = 0.04; BSinfo2(9,2) = 0.12; BSinfo2(9,3) = BSinfo2(9,2); 
BSinfo2(9,4) = 25; BSinfo2(9,5) = 0; BSname{9} = 'PS-S1-T1-1';
BSinfo2(10,1) = 0.04; BSinfo2(10,2) = 0.13; BSinfo2(10,3) = BSinfo2(10,2); 
BSinfo2(10,4) = 20; BSinfo2(10,5) = 0; BSname{10} = 'PS-S1-T1-2';
BSinfo2(11,1) = 0.055; BSinfo2(11,2) = 0.13; BSinfo2(11,3) = BSinfo2(11,2); 
BSinfo2(11,4) = 21; BSinfo2(11,5) = 0; BSname{11} = 'PS-S1-T2-1';
BSinfo2(12,1) = 0.05; BSinfo2(12,2) = 0.13; BSinfo2(12,3) = BSinfo2(12,2); 
BSinfo2(12,4) = 24; BSinfo2(12,5) = 0; BSname{12} = 'PS-S1-T2-2';

% Prestrained SMA,Araldite 2015
BSinfo2(13,1) = 0.08; BSinfo2(13,2) = 0.4; BSinfo2(13,3) = 0.6; 
BSinfo2(13,4) = 18; BSinfo2(13,5) = 18; BSname{13} = 'PS-A-T0.5-1';
BSinfo2(14,1) = 0.06; BSinfo2(14,2) = 0.3; BSinfo2(14,3) = 0.8; 
BSinfo2(14,4) = 17; BSinfo2(14,5) = 17; BSname{14} = 'PS-A-T0.5-2';
BSinfo2(15,1) = 0.07; BSinfo2(15,2) = 0.4; BSinfo2(15,3) = 0.8; 
BSinfo2(15,4) = 14; BSinfo2(15,5) = 14; BSname{15} = 'PS-A-T1-1';
BSinfo2(16,1) = 0.07; BSinfo2(16,2) = 0.6; BSinfo2(16,3) = 1; 
BSinfo2(16,4) = 13; BSinfo2(16,5) = 13; BSname{16} = 'PS-A-T1-2';
BSinfo2(17,1) = 0.1; BSinfo2(17,2) = 0.5; BSinfo2(17,3) = 0.8; 
BSinfo2(17,4) = 14; BSinfo2(17,5) = 14; BSname{17} = 'PS-A-T2-1';
BSinfo2(18,1) = 0.1; BSinfo2(18,2) = 0.4; BSinfo2(18,3) = 0.8; 
BSinfo2(18,4) = 16; BSinfo2(18,5) = 16; BSname{18} = 'PS-A-T2-2';

% Prestrained SMA,Sika 1277
BSinfo2(19,1) = 0.04; BSinfo2(19,2) = 0.4; BSinfo2(19,3) = 0.7; 
BSinfo2(19,4) = 20; BSinfo2(19,5) = 20; BSname{19} = 'PS-S2-T0.5-1';
BSinfo2(20,1) = 0.05; BSinfo2(20,2) = 0.5; BSinfo2(20,3) = 0.8; 
BSinfo2(20,4) = 20; BSinfo2(20,5) = 20; BSname{20} = 'PS-S2-T0.5-2';
BSinfo2(21,1) = 0.07; BSinfo2(21,2) = 0.4; BSinfo2(21,3) = 0.66; 
BSinfo2(21,4) = 20; BSinfo2(21,5) = 20; BSname{21} = 'PS-S2-T1-1';
BSinfo2(22,1) = 0.06; BSinfo2(22,2) = 0.6; BSinfo2(22,3) = 1.1; 
BSinfo2(22,4) = 16; BSinfo2(22,5) = 16; BSname{22} = 'PS-S2-T1-2';
BSinfo2(23,1) = 0.12; BSinfo2(23,2) = 0.83; BSinfo2(23,3) = 1; 
BSinfo2(23,4) = 20; BSinfo2(23,5) = 20; BSname{23} = 'PS-S2-T2-1';
BSinfo2(24,1) = 0.12; BSinfo2(24,2) = 0.86; BSinfo2(24,3) = 1.4; 
BSinfo2(24,4) = 16; BSinfo2(24,5) = 16; BSname{24} = 'PS-S2-T2-2';
% BSinfo2(23,:) = BSinfo2(24,:); BSname{23} = BSname{24}; 
%'PS-S2-T2-1' did not have information, therefore, use 'PS-S2-T2-2' to represent

% Prestrained SMA,Sika 888
BSinfo2(25,1) = 0.05; BSinfo2(25,2) = 0.17; BSinfo2(25,3) = 1; 
BSinfo2(25,4) = 30; BSinfo2(25,5) = 4; BSname{25} = 'NS-S3-T0.5-1';
BSinfo2(26,1) = 0.07; BSinfo2(26,2) = 0.25; BSinfo2(26,3) = 0.7; 
BSinfo2(26,4) = 18; BSinfo2(26,5) = 18; BSname{26} = 'PS-S3-T0.5-1';

% Parametric study
BS_pm(1,1) = 0.05; BS_pm(1,2) = 0.8; BS_pm(1,3) = BS_pm(1,2); 
BS_pm(1,4) = 20; BS_pm(1,5) = 0; BS_pm_name{1} = 'i. Original triangle'; %case 1, original triangle
BS_pm(2,1) = 0.05; BS_pm(2,2) = 0.3; BS_pm(2,3) = 0.55; 
BS_pm(2,4) = 20; BS_pm(2,5) = 20; BS_pm_name{2} = 'ii. Trapezoid with the same area'; %case 2, trapezoid with the same area, stiffness, and shear strength
BS_pm(3,1) = 0.3; BS_pm(3,2) = 0.8; BS_pm(3,3) = BS_pm(3,2); 
BS_pm(3,4) = 20; BS_pm(3,5) = 0; BS_pm_name{3} = 'iii. Triangle with reduced stiffness'; %case 3, triangle with reduced stiffness
BS_pm(4,1) = 0.1; BS_pm(4,2) = 0.4; BS_pm(4,3) = BS_pm(4,2); 
BS_pm(4,4) = 40; BS_pm(4,5) = 0; BS_pm_name{4} = 'iv. Triangle with increased strength'; %case 4, triangle with increased strength
BS_pm(5,1) = 0.05; BS_pm(5,2) = 0.4; BS_pm(5,3) = BS_pm(5,2); 
BS_pm(5,4) = 20; BS_pm(5,5) = 0; BS_pm_name{5} = 'v. Triangle with reduced area'; %case 5, triangle with reduced area

end