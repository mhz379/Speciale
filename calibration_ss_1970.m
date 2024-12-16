function [obj_val] = calibration_ss_1970(param_vector,moments)

global oo_ D_par_1970 theta_par_1970 alpha_par CDY

D_par_1970 = param_vector(1);
filename='data/D_1970.xlsx';
writematrix(D_par_1970,filename);
theta_par_1970 = param_vector(2);
filename='data/theta_1970.xlsx';
writematrix(theta_par_1970,filename);
alpha_par = param_vector(3);
filename='data/alpha.xlsx';
writematrix(alpha_par,filename);
beta_par = param_vector(4);
filename='data/beta.xlsx';
writematrix(beta_par,filename);


%try
    dynare dynare_ss_1970_calibr nostrict
%catch
%    disp('errore')
%end

iter=1;

% Generate the objective function: this objective function will be minimized during the calibration of the 1970 economy 
obj_val = 50000*(oo_.steady_state(344,1) - moments.r)^2 + 1000*(oo_.steady_state(359,1) - moments.ls)^2 + 1000*(CDY - moments.debt_inc)^2 + 1000*(oo_.steady_state(358,1) - moments.IY)^2;


obj_val

end