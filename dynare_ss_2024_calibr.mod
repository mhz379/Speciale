
%--------------------------------------------------------------------------
% Authors: Alex Crescentini & Federico Giri.
% Università Politecnica delle Marche, Ancona, Italy.
% January 2023.
%--------------------------------------------------------------------------
% Dynare Replication Code for:
% "A Model of Secular Stagnation: Theory and Quantitative Evaluation"
% Eggertsson, Mehrotra and Robbins (2019), American Economic Journal.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
//DESCRIPTION
%--------------------------------------------------------------------------

%Code used for calibration targets at 2024

%--------------------------------------------------------------------------
//LOAD DATA TO BE USED IN THE MODEL
%--------------------------------------------------------------------------

ss_2024 = xlsread('data/ss_2024_calibr.xlsx');

%Parameters to be Calibrated at 2024 Targets
D_par_par = xlsread('data/D_2024.xlsx');
theta_par_par = xlsread('data/theta_2024.xlsx');
mu_par_par = xlsread('data/mu.xlsx');

%Parameters to be Calibrated at 1970 Targets
alpha_par_par = xlsread('data/alpha.xlsx');
beta_par_par = xlsread('data/beta.xlsx');

%--------------------------------------------------------------------------
//GLOBAL VARIABLES
%--------------------------------------------------------------------------

@#define J = 74
@#define T = 201

%--------------------------------------------------------------------------
//ENDOGENOUS VARIABLES
%--------------------------------------------------------------------------

var 

@#for j in 1:J
	n@{j}
@#endfor

@#for j in 1:J
	c@{j}
@#endfor

@#for j in 1:J+1
	a@{j}
@#endfor

@#for j in 1:J
	lambda@{j}
@#endfor

@#for j in 1:J-32
	pi@{j}
@#endfor

x74, q25, w, rk, r, price, tau, PI, Y, N, L, C, K, gov_rev, gov_deficit
gov_debt, I, A_adj, IY, LS, BY, E, d 

;

%--------------------------------------------------------------------------
//PARAMETERS
%--------------------------------------------------------------------------

parameters

alpha, delta, beta, gamma, sigma, mu, g, d_bar, taup,

@#for j in 1:J-32
	hc@{j}
@#endfor

@#for j in 1:J
    s@{j}
@#endfor

@#for j in 1:J
    sv@{j}
@#endfor

@#for j in 1:J
    su@{j}
@#endfor

n, fert_25, e, b, D, theta, AL_growth, AL, AK
    
;

%--------------------------------------------------------------------------
//PARAMETERS' CALIBRATION
%-------------------------------------------------------------------------- 

%%Parameters that do not change with SS
	@#for j in 1:42
        hc@{j}=datahc(1,@{j});
    @#endfor
    delta=delta_p;
    gamma=gamma_p;
    sigma=sigma_p; 
    g=g_p;
    d_bar = d_bar_p;
    taup = taup_p;

%%Parameters that change with SS
    @#for j in 1:J
        s@{j}=datasurvival1_2024(1,@{j});
    @#endfor
    @#for j in 1:J
        sv@{j}=datasurvival2_2024(1,@{j});
    @#endfor
    @#for j in 1:J
        su@{j}=data_uncond_survival_2024(1,@{j});
    @#endfor     
    n=n_2024;               
    fert_25=fert_25_2024;                                    
    b=b_2024;               
    AL_growth=AL_growth_2024;
    e=1;       
    AL=1;                             
    AK=1;  
                               

%To be Calibrated at 2024 Targets
    mu=mu_par_par;
    D=D_par_par;
    theta=theta_par_par;

%To be Calibrated at 1970 Targets
    alpha=alpha_par_par;
    beta=beta_par_par;

%--------------------------------------------------------------------------
//MODEL BLOCK: EQUATIONS
%--------------------------------------------------------------------------

model;

//Demographics
%--------------------------------------------------------------------------

[name = 'Demographics Equation 1']
n1=1;

@#for j in 1:J-1
    [name = 'Demographics Equation @{j+1}']
    n@{j+1} = (s@{j}*n@{j})/(1+n);
@#endfor

//Households: First-Order conditions
%--------------------------------------------------------------------------

@#for j in 1:J-1
    [name = 'Households First-order Condition @{j}']
    (1/beta) = ((c@{j+1}/c@{j})^(-1/gamma))*(1+r) + (lambda@{j+1})*((c@{j})^(1/gamma))/(su@{j}*beta^@{j}*e); %har fjernet sv@{j+1} i tælleren.
@#endfor

    [name = 'Households First-order Condition 74']
    x74 = (fert_25^2/mu)^(-gamma)*c74; %opløftet fert_25 i anden pga. børnebørn.
    
//Households: Budget Constraints
%--------------------------------------------------------------------------

[name = 'Households Budget Constraint: initial condition']
a1=0;

@#for j in 1:J-51
    [name = 'Households Budget Constraint @{j}']
    a@{j+1} = (1+r)/(sv@{j})*a@{j} + ((1-taup)*pi@{j}+(1-tau)*w*hc@{j})*(1+AL_growth)^@{j} - c@{j};
@#endfor

[name = 'Households Budget Constraint 24']
a25 = (1+r)/sv24*a24 + ((1-taup)*pi24+(1-tau)*w*hc24)*(1+AL_growth)^@{24} + q25*(1+AL_growth)^@{25} - c24;

[name = 'Bequest equation']
q25 = ((x74*fert_25^2*n74)/(n25)*(1/(1+AL_growth)^75))*1/(1+n); %opløftet fert_25 i anden pga. børnebørn.

@#for j in 25:J-32
    [name = 'Households Budget Constraint @{j}']
    a@{j+1} = (1+r)/sv@{j}*a@{j} + ((1-taup)*pi@{j}+(1-tau)*w*hc@{j})*(1+AL_growth)^@{j} - c@{j};
@#endfor

@#for j in 43:J-1
    [name = 'Households Budget Constraint @{j}']
    a@{j+1} = (1+r)/sv@{j}*a@{j} - c@{j} + d*w*(1+AL_growth)^@{j}; %skalerer dw
@#endfor
    
    [name = 'Households Budget Constraint 74']
    a75 = (1+r)/sv74*a74 - (fert_25^2*x74) - c74 + d*w*(1+AL_growth)^@{74}; %opløftet fert_25 i anden pga. børnebørn. Skalerer dw.

[name = 'Households Budget Constraint: terminal condition']
a75 = 0;

//Households: Financial Constraints (OBCs)
%--------------------------------------------------------------------------

@#for j in 1:J-32
    [name = 'Households Financial Constraint @{j}']
    min(lambda@{j},a@{j}+(D*w*hc@{j})*(1+AL_growth)^@{j}) = 0;  
@#endfor
@#for j in 43:J
    [name = 'Households Financial Constraint @{j}']
    min(lambda@{j},a@{j}) = 0;  
@#endfor

//Households: Profits from Firms
%--------------------------------------------------------------------------

@#for j in 1:J-32
    [name = 'Households Profits from Firms @{j}']
    pi@{j}=hc@{j}*PI/L;
@#endfor

//Firms
%--------------------------------------------------------------------------

price = (theta-1)/theta;
A_adj = (price*(((alpha)*(AK*K)^((sigma-1)/sigma)+(1-alpha)*(AL*L)^((sigma-1)/sigma))^(1/(sigma-1)))*(1-alpha)*AL^((sigma-1)/sigma)*L^(-1/sigma));
w  = 1; %Er lig 1, da AL = 1.
rk = (price*(((alpha)*(AK*K)^((sigma-1)/sigma)+(1-alpha)*(AL*L)^((sigma-1)/sigma))^(1/(sigma-1)))*(alpha)*AK^((sigma-1)/sigma)*K^(-1/sigma))/A_adj;
Y  = (((alpha)*(AK*K)^((sigma-1)/sigma)+(1-alpha)*(AL*L)^((sigma-1)/sigma))^(sigma/(sigma-1)))/A_adj;
r = rk/e + (1-delta)*e/e - 1;
PI = Y/theta;

 //Government
%--------------------------------------------------------------------------

[name = 'Government']
E = (
@#for j in 43:J
    + n@{j}*d*w
@#endfor
);

d = d_bar*(1-tau); %kompensationsgrad før skat

gov_deficit*gov_rev = ((1+AL_growth)*(1+n)-1)*(gov_debt*K);
gov_debt=b*Y/K;
gov_rev = (g*Y+E+r*gov_debt*K);
tau*w = ((gov_rev)*(1-gov_deficit))/(L) - taup*PI/L;

//Aggregates
%--------------------------------------------------------------------------

[name = '(N) Population']
N = (
@#for j in 1:J
    + n@{j}
@#endfor
);

[name = '(L) Labor']
L = (
@#for j in 1:J-32
    + n@{j}*hc@{j}
@#endfor
);

[name = '(C) Consumption']
C = (
@#for j in 1:J
    + n@{j}*c@{j}/(1+AL_growth)^(@{j})
@#endfor
);

[name = '(K) Kapital']
K = ((
@#for j in 1:J
    + n@{j}*a@{j}/sv@{j}/(1+AL_growth)^(@{j}) %dividerer med sv@{j} da kapitalen er prædetermineret
@#endfor
)/(e+gov_debt));

[name = '(I) Investment']
I = (1+AL_growth)*(1+n)*e*K - (1-delta)*e*K;

[name = 'Investment Output Ratio (I/Y)']
IY = I/Y;

[name = 'Labor Share']
LS = L*w/(L*w+r*K+PI+delta*K); %har ganget w på L (selvom w=1)

[name = 'Bequest to output']
BY=q25*n25/Y;

end;

%--------------------------------------------------------------------------
//STEADY-STATE BLOCK: Initial Values
%--------------------------------------------------------------------------

initval;

@#for j in 1:J
    n@{j}=ss_2024(@{j});
@#endfor
    
@#for j in 1:J
	c@{j}=ss_2024(@{j+74});
@#endfor
        
@#for j in 1:J+1
	a@{j}=ss_2024(@{j+148});
@#endfor
            
@#for j in 1:J
	lambda@{j}=ss_2024(@{j+223});
@#endfor
                
@#for j in 1:J-32
	pi@{j}=ss_2024(@{j+297}); 
@#endfor

%disse skal også ændres hvis J ændres
x74=ss_2024(340); 
q25=ss_2024(341); 
w=ss_2024(342); 
rk=ss_2024(343); 
r=ss_2024(344); 
price=ss_2024(345); 
tau=ss_2024(346); 
PI=ss_2024(347); 
Y=ss_2024(348); 
N=ss_2024(349); 
L=ss_2024(350); 
C=ss_2024(351);
K=ss_2024(352); 
gov_rev=ss_2024(353); 
gov_deficit=ss_2024(354); 
gov_debt=ss_2024(355); 
I=ss_2024(356);
A_adj=ss_2024(357);
IY=ss_2024(358);
LS=ss_2024(359);
BY=ss_2024(360); 
E = ss_2024(361);
d = ss_2024(362);

end;

%--------------------------------------------------------------------------
//STEADY-STATE BLOCK: Computation
%--------------------------------------------------------------------------

options_.debug=1;
%options_.dynatol.f=1e-3;
steady(maxit=1000, solve_algo=0);
%resid;
%check;

%--------------------------------------------------------------------------
//NEEDED FOR CALIBRATION
%--------------------------------------------------------------------------

% Load the AL_growth_2024 value from the .mat file
load('AL_growth_2024_value.mat', 'AL_growth_2024');
%AL_growth= 0.0085;

capi=oo_.steady_state(149:223,1); %ændret fra 113:169:1. Skal muligvis være 224 i stedet for 223.
filename='data/ss_2024_calibr.xlsx';
writematrix(oo_.steady_state,filename,'WriteMode','replacefile');

%Consumer Debt to GDP
datasurvival2_2024 = readmatrix('data/data_survival_2_2024.xlsx');
popgen=oo_.steady_state(1:74,1);
output=oo_.steady_state(348,1);
personal_debt=0;
J=74;
for j = 1:J
	a_scaled2_2024(j,1)=capi(j,1)/(1+AL_growth)^(j);
end
for j = 1:J
	personal_debt=personal_debt+(a_scaled2_2024(j,1)<0)*a_scaled2_2024(j,1)*popgen(j,1)/datasurvival2_2024(1,j);
end
global CDY
CDY=-personal_debt/output;
filename='data/CDY_2024.xlsx';
writematrix(CDY,filename,'WriteMode','replacefile');

%--------------------------------------------------------------------------
//LATEX FILES
%--------------------------------------------------------------------------

%write_latex_dynamic_model;
%write_latex_static_model;
%write_latex_definitions;
%write_latex_parameter_table;
%collect_latex_files;