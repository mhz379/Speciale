%--------------------------------------------------------------------------
% VALDER Ø. FREDENS, SPECIALE / MASTER'S THESIS
% Kode til simularing af baseline / code to simulate the baseline
% 
% Koden er baseret på / the code is based on:
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

%Code used for calibration targets at 1970

%--------------------------------------------------------------------------
//LOAD DATA TO BE USED IN THE MODEL
%--------------------------------------------------------------------------

ss_1970 = xlsread('data/ss_1970.xlsx');

%Calibrated Parameters at 2024 Targets
calibrated_parameters_2024=xlsread('calibration_2024_baseline_dynare.xlsx');
%Parameters to be Calibrated at 1970 Targets
D_par_par=xlsread('data/D_1970.xlsx');
theta_par_par=xlsread('data/theta_1970.xlsx');
alpha_par_par=xlsread('data/alpha.xlsx');
beta_par_par=xlsread('data/beta.xlsx');

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

x74, q25, w, rk, r, price, tau, PI, Y, N, L, C, K, gov_rev, gov_deficit,
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
        s@{j}=datasurvival1_1970(1,@{j});
    @#endfor
    @#for j in 1:J
        sv@{j}=datasurvival2_1970(1,@{j});
    @#endfor
    @#for j in 1:J
        su@{j}=data_uncond_survival_1970(1,@{j});
    @#endfor 
    n=n_1970;               
    fert_25=fert_25_1970;                                        
    b=b_1970;                                   
    AL_growth=AL_growth_1970; 
    e=1.0;  
    %e = load('e_1970_value.mat', 'e_1970');       
    AL=1;                             
    AK=1;                             

%%Parameters calibrated for 2024 Targets
    mu=calibrated_parameters_2024(3);
   
%%Parameters to be calibrated for 1970 Targets
	D=D_par_par; 
	theta=theta_par_par;
    alpha = alpha_par_par;
    beta = beta_par_par;
    
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

[name = 'Households Budget Constraint: Initial Condition']
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
    a@{j+1} = (1+r)/sv@{j}*a@{j} - c@{j} + d*w*(1+AL_growth)^@{j}; %skalerer dw;
@#endfor
    
[name = 'Households Budget Constraint 74']
a75 = (1+r)/sv74*a74 - (fert_25^2*x74) - c74 + d*w*(1+AL_growth)^@{74}; %fertilitet opløftet i anden pga. børnebørn. %skalerer dw
 
[name = 'Households Budget Constraint: Terminal Condition']
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
w  = 1;
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
LS = L/(L+r*K+PI+delta*K);

[name = 'Bequest to output']
BY=q25*n25/Y;

end;

%--------------------------------------------------------------------------
//STEADY-STATE BLOCK: Initial Values
%--------------------------------------------------------------------------

initval;

@#for j in 1:J
    n@{j}=ss_1970(@{j});
@#endfor

@#for j in 1:J
	c@{j}=ss_1970(@{j+74}); 
@#endfor

@#for j in 1:J+1
	a@{j}=ss_1970(@{j+148});
@#endfor

@#for j in 1:J
	lambda@{j}=ss_1970(@{j+223});
@#endfor

@#for j in 1:J-32 
	pi@{j}=ss_1970(@{j+297});
@#endfor

%disse skal også ændres, hvis J ændres.
x74=ss_1970(340);
q25=ss_1970(341);
w=ss_1970(342);
rk=ss_1970(343);
r=ss_1970(344);
price=ss_1970(345);
tau=ss_1970(346);
PI=ss_1970(347); %herfra og ned passer ikke med ss_1970
Y=ss_1970(348); 
N=ss_1970(349); 
L=ss_1970(350);
C=ss_1970(351);
K=ss_1970(352);
gov_rev=ss_1970(353);
gov_deficit=ss_1970(354);
gov_debt=ss_1970(355);
I=ss_1970(356);
A_adj=ss_1970(357);
IY = ss_1970(358); %tilføjet
LS = ss_1970(359); %tilføjet
BY = ss_1970(360);
E=ss_1970(361); %tilføjet
d=ss_1970(362); %tilføjet

end;

%--------------------------------------------------------------------------
//STEADY-STATE BLOCK: Computation
%--------------------------------------------------------------------------

%options_.debug=1;
steady(maxit=1000, solve_algo=0);
%resid;
%check;

%--------------------------------------------------------------------------
//OUTPUT AS IN AUERBACH & KOTLIKOFF (1987) AND EGGERTSSON ET AL. (2019)
%--------------------------------------------------------------------------

%disp(oo_.steady_state(340:352,1)); %behøves ikke

disp(['   '])
disp(['DYNARE-Statistics'])
disp(['   '])
disp(['Year  ' 'Capital (K)  ' 'Labor (L)  ' 'Population (N)'])
disp([num2str(0) '     ' num2str(oo_.steady_state(352,1)) '       ' num2str(oo_.steady_state(350,1)) '     ' num2str(oo_.steady_state(349,1))])

disp(['   '])
disp (['Income (Y)  ' 'Consumption (C)  ' 'Agg. Profit (PI)'])
disp([num2str(oo_.steady_state(348,1)) '     ' num2str(oo_.steady_state(351,1)) '          ' num2str(oo_.steady_state(347,1)) '    '])

disp(['   '])
disp(['Wage (w)  ' 'Rental K (rk)  ' 'Interest (r)  '   'Wage Tax (tau)'])
disp([num2str(oo_.steady_state(342,1)) '         ' num2str(oo_.steady_state(343,1)) '        ' num2str(oo_.steady_state(344,1)) '      ' num2str(oo_.steady_state(346,1))]) %ændret

disp(['   '])
disp(['Bequest (q25*fert_25)  ' 'q25     ' 'x74/(1+AL_growth)^74'])
disp([num2str(oo_.steady_state(341,1)*fert_25) '                ' num2str(oo_.steady_state(341,1)) '  ' num2str(oo_.steady_state(341,1)/(1+AL_growth)^74)]) %ændret, også potenserne

disp(['   '])
disp(['Pop Growth (n)  '])
disp([num2str(n)])

disp(['   '])
disp(['Debt (b*Y/K)  '])
disp([num2str(b*(oo_.steady_state(348,1)/oo_.steady_state(352,1)))]) %ændret

disp(['   '])
disp(['G (g*Y)  '])
disp([num2str(g*oo_.steady_state(348,1))]) %ændret

%--------------------------------------------------------------------------
//FILES PREPERATION FOR TRANSITION DYNAMICS: "dynare_transition.mod"
%--------------------------------------------------------------------------

%%1970: SCALING FOR PRODUCTIVITY GROWTH%%

% Load the AL_growth_2024 value from the .mat file
load('AL_growth_1970_value.mat', 'AL_growth_1970');

J=74;
%AL_growth=0.020209202089143; 

cons=oo_.steady_state(75:148,1); %ændret. Skal muligvis være 75
capi=oo_.steady_state(149:223,1); %ændret fra 113:169:1
for j = 1:J
	c_scaled_1970(j,1)=cons(j,1)/(1+AL_growth)^(j);
end
filename='data/c_scaled_1970.xlsx';
writematrix(c_scaled_1970,filename,'WriteMode','replacefile');
J=75;
for j = 1:J
	a_scaled_1970(j,1)=capi(j,1)/(1+AL_growth)^(j-1);
end
filename='data/a_scaled_1970.xlsx';
writematrix(a_scaled_1970,filename,'WriteMode','replacefile');
x74_scaled_1970=oo_.steady_state(340,1)/(1+AL_growth)^74;
filename='data/x74_scaled_1970.xlsx';
writematrix(x74_scaled_1970,filename,'WriteMode','replacefile');
filename='data/ss_1970.xlsx';
writematrix(oo_.steady_state,filename,'WriteMode','replacefile');

%Consumer Debt to GDP
datasurvival2_1970 = readmatrix('data/data_survival_2_1970.xlsx');
popgen=oo_.steady_state(1:74,1);
output=oo_.steady_state(348,1);
personal_debt=0;
J=74;
for j = 1:J
	a_scaled2_1970(j,1)=capi(j,1)/(1+AL_growth)^(j);
end
for j = 1:J
	personal_debt=personal_debt+(a_scaled2_1970(j,1)<0)*a_scaled2_1970(j,1)*popgen(j,1)/datasurvival2_1970(1,j);
end
global CDY
CDY=-personal_debt/output;
filename='data/CDY_1970.xlsx';
writematrix(CDY,filename,'WriteMode','replacefile');

%--------------------------------------------------------------------------
//LATEX FILES
%--------------------------------------------------------------------------

%write_latex_dynamic_model;
%write_latex_static_model;
%write_latex_definitions;
%write_latex_parameter_table;
%collect_latex_files;