function [obj_val] = calibration_ss_2024(param_vector,moments)

global oo_ D_par_2024 theta_par_2024 mu_par CDY  

%alpha_par = param_vector(1);
%filename='data/alpha.xlsx';
%writematrix(alpha_par,filename);
%beta_par = param_vector(1); 
%filename='data/beta.xlsx';
%writematrix(beta_par,filename);


D_par_2024 = param_vector(1);
filename='data/D_2024.xlsx';
writematrix(D_par_2024,filename);
theta_par_2024 = param_vector(2);
filename='data/theta_2024.xlsx';
writematrix(theta_par_2024,filename);
mu_par = param_vector(3);
filename='data/mu.xlsx';
writematrix(mu_par,filename);

moments_2024.beq_inc = .005;

%%Sometimes try and catch is useful in making the calibration with Dynare
%try
    dynare dynare_ss_2024_calibr
%catch
 %   disp('SS not found trying a new guess without stopping the algorithm')
%end

iter=1;

% Generate the objective function: this objective function will be minimized during the calibration of the 2024 economy 
obj_val = 1000*(oo_.steady_state(359,1) - moments.ls)^2 + 1000*(CDY - moments.debt_inc)^2 + 1000*(oo_.steady_state(360,1) - moments_2024.beq_inc)^2;

obj_val

end