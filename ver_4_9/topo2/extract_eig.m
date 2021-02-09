clc 
clear 
close all 

% load('ch64_win64_seg64_frL25_frU25_overLap25.mat')
load('ch64_win64_seg64_frL12_frU12_overLap25.mat')

eig_val = sorted_eig_info{1};
eig_vec = sorted_eig_info{2};

Eig_Info.eig_val = eig_val;
Eig_Info.eig_vec = eig_vec;

save('Eig_Info_12.mat', 'Eig_Info');
k