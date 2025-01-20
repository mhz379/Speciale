%--------------------------------------------------------------------------
% VALDER Ø. FREDENS, SPECIALE / MASTER'S THESIS
% Kode til simularing af baseline / code to simulate the baseline
% % Koden er baseret på / the code is based on:
% Authors: Alex Crescentini & Federico Giri.
% Università Politecnica delle Marche, Ancona, Italy.
% January 2023.
%--------------------------------------------------------------------------
% Dynare Replication Code for:
% "A Model of Secular Stagnation: Theory and Quantitative Evaluation"
% Eggertsson, Mehrotra and Robbins (2019), American Economic Journal.
%--------------------------------------------------------------------------

%% DESCRIPTION

%The following file compares the output of the model coming from the
%original Matlab code by Eggertsson, Mehrotra and Robbins (2019) and
%by the replication with Dynare from us.

%In detail: First it runs the main calibration and produces the 
%steady state of the model at 1970, at 2024 and the transitional dynamics,
%with both our Dynare code (point A) and the original Matlab code (point B) 
%of the authors. Second (point C), it compares the output with the one 
%obtained by Eggertsson, Mehrotra and Robbins (2019).

%In doing so, you need to change your path and choose if you want to
%recalibrate the model or not.

%% A) DYNARE (our replication code)

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------BASELINE CALIBRATION--------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_calibration=0; %Choose run_calibration=1 if you want to recalibrate the model, otherwise load the already calibrated parameters

%%Parameters that do not change with SS

    datahc = readmatrix('data/data_hc.xlsx');     %Humankapitalprofil. Fra NTA
    delta_p=0.10;                                 %Egg: 0.1244. 
    gamma_p=0.95;                                 %Intertemporal sub. elasticitet for forbrug. Hedder rho i artiklen. Egg: 0.75. 
    sigma_p=0.6;                                  %Sub. elasticitet for K/L: Hvis sigma = 1, bliver den Cobb-Douglas, og så skal der ændres nogle ligninger
    g_p = 0.1554;                                 %Off. forbrug som andel af BNP  ekskl. pensioner og investeringer og renteudgifter. Gns. fra NIPA. Egg.: 0.2128
    d_bar_p = 0.505; %                            %Kompensationsgrad efter skat: 0,505. Fra OECD. 
    taup_p = 0.15;                                %Skat på profit. Empirisk eff. capital gains tax på 15%.

%%Parameters that change with SS

    %2024 (det vil sige de parametre, som skal gælde, når modellen slutter i steady state i den terminale periode)
    datasurvival1_2024 = readmatrix('data/data_survival_1_2024.xlsx');                  %2024 overlevelsessandsynligheder. Fra FN. 
    datasurvival2_2024 = readmatrix('data/data_survival_2_2024.xlsx');                  %2024 overlevelsessandsynligheder. 
    data_uncond_survival_2024 = readmatrix('data/data_uncondit_survival_2024.xlsx');    %2024 overlevelsessandsynligheder.
    n_2024 = -0.007906615;                                                              %Befolkningsvækst i steady state. Baseret på ligning (2). Beregnet i exo_matrix_start100. Egg: -0.002578211534468.
    fert_25_2024 = 0.82;                                                                %2024 Fertilitet i ss (fra 2021). Se FN og exo_matrix. Egg: 0.9375.
    b_2024 = 1.2;                                                                       %2024 statsgæld i % af BNP i 2021 (gns.) (se rådata).
    e_2024 = 1;                                                                         %2024 relativ pris på kapital. Bruges ikke. Egg.: 1.009368277002316.
    AL_growth_2024 = 0.0085;                                                            %2024 produktivitetsvækst (gns. inden 2021. Se Ny data). Egg: 0.006460448230438.
    AL_2024 = 5.465873263;                                                              %Slutværdi for produktivitetsfaktor (niveau). Sidder i celle E201 i exo_matrix_start100. Egg: 2.99049071248668.                                                            
    AK_2024=1;                                                                          %Antages at være lig 1 hele tiden, jf. Harrod-neutralitet.
        
    %1970
    datasurvival1_1970 = readmatrix('data/data_survival_1_1970.xlsx');                  %1970 overlevelsessandsynligheder. Fra FN.
    datasurvival2_1970 = readmatrix('data/data_survival_2_1970.xlsx');                  %1970 overlevelsessandsynligheder 
    data_uncond_survival_1970 = readmatrix('data/data_uncondit_survival_1970.xlsx');    %1970 overlevelsessandsynligheder
    n_1970= 0.013549868016229;                                                          %1970 befolkningsvækst i 1970 ss (se fertilitetsdata).
    fert_25_1970=1.4;                                                                   %1970 fertilitetskvotient (se fertilitetsdata).
    b_1970= 0.35;                                                                       %1970 statsgæld i % af BNP (se rådata). Egg.: 0.403662983955455.
    e_1970= 1.0;                                                                        %1970 relativ pris på kapital. Tager den ikke med. Egg.: 1.3.
    AL_growth_1970= 0.02;                                                               %1970 produktivitetsvækst i ss (see rådata). %Egg.: 0.020209202089143.
    AL_1970=1;                                                                          %Normaliseres til at starte på 1.
    AK_1970=1;                                                                          %Er hele tiden lig 1 jf. Harrod-neutralitet.
    
% Save the value of AL_growth_2024 into a file or workspace that can be accessed by dynare_ss_2024
save('AL_growth_2024_value.mat', 'AL_growth_2024');
save('AL_growth_1970_value.mat', 'AL_growth_1970');
save('e_1970_value.mat', 'e_1970'); %Hvis den relative pris på kapital er med.

 
if run_calibration==1

   %%CALIBRATION for 1970 TARGETS
    
        %Targets: these are the moments we match for 1970
        moments_1970.debt_inc = .02240;  %Egg: 0.0421. Drevet af gældskvoten D.
        moments_1970.ls = 0.628;         %Egg.: 0.7240. 62,8% er det glidende gns. fra BLS.
        moments_1970.IY = 0.2167;        %Drevet af alpha. I/Y=0.217 i 1970 (glidende gns.)
        moments_1970.r = 0.043;          %Drevet af beta. 4,3% er den udglattede værdi fra HLW.
        %Parameters bounds and guesses
        lb = [.05;2; 0.25;0.96];
        ub = [.7;20; 0.4;1];
        param_guess = [0.2;6;0.37;0.98];
        %Run Calibration
        options = optimoptions('fmincon','Display','iter','DiffMinChange',0.00005,'OptimalityTolerance',1e-10);
        [calibration_1970_baseline_dynare] = fmincon(@(param) calibration_ss_1970(param,moments_1970),param_guess,[],[],[],[],lb,ub,[],options);
        %[calibration_1970_baseline_dynare] = fminsearch(@(param) calibration_ss_1970(param,moments_1970),param_guess);
        %Save Results
        writematrix(calibration_1970_baseline_dynare,'calibration_1970_baseline_dynare.xlsx');
        calibrated_parameters_1970=xlsread('calibration_1970_baseline_dynare.xlsx');
    
    %%CALIBRATION for 2024 TARGETS.
    
        %Targets: these are the moments we match for 2024. 
        %moments_2024.IY = .1822; %IY i 2021. 0.159 før; %driven by alpha
        %moments_2024.r = 0.00; %ændret fra -.0147;  %driven by beta
        moments_2024.debt_inc = 0.12345;    %Drevet af D. Egg: 0.0633. Glidende gns. for 2015 (se noter).  
        moments_2024.ls = 0.562;         %Drevet af theta. Glidende gns. fra BLS. Egg.: 0.6599.
        moments_2024.beq_inc = .005;     %Drevet af mu. Egg.: 0.03. Vælger 0,5% manuelt.
        %Parameters bounds and guesses
        lb = [.05;3;10];
        ub = [.99;12;13.6]; 
        param_guess = [0.7;3;13.4]; 
        %Run Calibration
        [calibration_2024_baseline_dynare] = fmincon(@(param) calibration_ss_2024(param,moments_2024),param_guess,[],[],[],[],lb,ub);
        %[calibration_2024_baseline_dynare] = fminsearch(@(param) calibration_ss_2024(param,moments_2024),param_guess);
        %Save Results
        writematrix(calibration_2024_baseline_dynare,'calibration_2024_baseline_dynare.xlsx');
        calibrated_parameters_2024=xlsread('calibration_2024_baseline_dynare.xlsx');
        
        
        
else 
    
    calibrated_parameters_1970=xlsread('calibration_1970_baseline_dynare.xlsx');
    calibrated_parameters_2024=xlsread('calibration_2024_baseline_dynare.xlsx');

end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------SIMULATION---------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Steady State 1970
dynare dynare_ss_1970;
%Steady State 2024
dynare dynare_ss_2024;
%Transitional Dynamics
dynare dynare_transition nostrict;


%% B) MATLAB (original code from Eggertsson, Mehrotra and Robbins (2019))

clear all
close all

%%%%%%%%%%%%PUT YOUR PATH HERE AND BELOW%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We have set TOLERANCE (run_schedule.tol) = 1e-5 and MAXITER (run_schedule.maxiter) = 200
%Windows (Federico)
%run('D:\Programmer\Dynare\6.0\matlab\kode\Speciale\Eggertsson2019_CrescentiniGiri_2023\Dynare_Repl_Code\Eggertsson_Matlab_Code\114159-V1\data\Section-8\sec_stag_calibration_control_panel.m');
%Slå linje 141 fra, så den ikke genberegner Eggertsson.
%macOS (Alex)
%run('/....add your path of this folder here...../EI_ReplicEggertsson2019_CrescentiniGiri_2023/Dynare_Repl_Code/Eggertsson_Matlab_Code/114159-V1/data/Section-8/sec_stag_calibration_control_panel.m');


%% C) LOAD OUTPUT FROM DYNARE AND MATLAB CODES and PRODUCES PLOTS FOR COMPARISON
 
clear all
close all

%dynare output
load dynare_transition.mat 

%%%%%%%%%%%%PUT YOUR PATH HERE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%matlab output
%Windows
%load '.....add your path of this folder here......\EI_ReplicEggertsson2019_CrescentiniGiri_2023\Dynare_Repl_Code\Eggertsson_Matlab_Code\114159-V1\data\Section-8\matlab_transition.mat'
load 'D:\Programmer\Dynare\6.0\matlab\kode\Speciale\Eggertsson2019_CrescentiniGiri_2023\Dynare_Repl_Code\Eggertsson_Matlab_Code\114159-V1\data\Section-8\matlab_transition.mat'


%%Output Comparison
%Sammenligner resultaterne med Eggertsson m.fl.
run('plot_comparison.m');