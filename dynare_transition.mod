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

%Code for Transitional Dynamics

%--------------------------------------------------------------------------
//LOAD DATA TO BE USED IN THE MODEL
%--------------------------------------------------------------------------

ss_1970 = readmatrix('data/ss_1970.xlsx');
ss_2024 = readmatrix('data/ss_2024.xlsx'); 
exo_matrix = readmatrix('data/exo_matrix_start100.xlsx');
c_scaled_1970 = readmatrix('data/c_scaled_1970.xlsx');
a_scaled_1970 = readmatrix('data/a_scaled_1970.xlsx');
x74_scaled_1970 = readmatrix('data/x74_scaled_1970.xlsx');
c_scaled_2024 = readmatrix('data/c_scaled_2024.xlsx'); 
a_scaled_2024 = readmatrix('data/a_scaled_2024.xlsx'); 
x74_scaled_2024 = readmatrix('data/x74_scaled_2024.xlsx'); 

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

x74, q25, w, rk, r, price, tau, PI, Y, N, L, C, K, I, gov_rev, gov_deficit,
gov_debt, E, d

@#for j in 1:J-32
	opt_inc@{j}
@#endfor

;

%--------------------------------------------------------------------------
//EXOGENOUS VARIABLES
%--------------------------------------------------------------------------

varexo 

fert_25, e, b, theta, AL,

@#for j in 1:J
    s@{j}
@#endfor

@#for j in 1:J
    sv@{j}
@#endfor

@#for j in 1:J
    su@{j}
@#endfor

@#for j in 1:42
    D@{j}
@#endfor

;

%--------------------------------------------------------------------------
//PARAMETERS
%--------------------------------------------------------------------------

parameters

alpha, beta, delta, gamma, sigma, mu, d_bar, g, taup,

@#for j in 1:J-32
	hc@{j}
@#endfor

AK, A, A_adj

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
    d_bar=d_bar_p;
    AK=AK_1970; 
    A=1;
    A_adj = ss_1970(357); % ændret, IT COMES FROM THE 1970 SS
    taup = taup_p;
                                                                  
%%Parameters calibrated for 2024 Targets
    %alpha=calibrated_parameters_2024(1);
    %beta=calibrated_parameters_2024(1);
    mu=calibrated_parameters_2024(3);

%%Parameters calibrated for 1970 Targets
    alpha=calibrated_parameters_1970(3);
    beta=calibrated_parameters_1970(4);

%--------------------------------------------------------------------------
//MODEL BLOCK: EQUATIONS
%--------------------------------------------------------------------------

model;

//Demographics
%--------------------------------------------------------------------------

[name = 'Demographics Equation 1']
n1 = 1/su25(-1) * n25(-1) * fert_25; %n1: 25-årige i dag. n25(-1): forældrene i perioden før. Er 50 år i dag.

@#for j in 1:J-1
    [name = 'Demographics Equation @{j+1}']
    n@{j+1} = (s@{j}(-1)*n@{j}(-1));
@#endfor

//Households: First-Order conditions
%--------------------------------------------------------------------------

@#for j in 1:J-1
    [name = 'Households First-order Condition @{j}']
    (1/beta) = ((c@{j+1}(+1)/c@{j})^(-1/gamma))*(1+r(+1)) + (lambda@{j+1})*((c@{j})^(1/gamma))/(su@{j}*beta^@{j}*e); %har fjernet sv@{j+1} i tælleren.
@#endfor

    [name = 'Households First-order Condition 74']
    x74 = ((fert_25(-73)*fert_25(-48)/mu)^(-gamma))*c74; %fertilitet for både bedsteforældre selv og deres børns fertilitet.

//Households: Budget Constraints
%--------------------------------------------------------------------------

[name = 'Households Budget Constraint: initial condition']
a1=0;

@#for j in 1:J-51
    [name = 'Households Budget Constraint @{j}']
    a@{j+1} = (1+r)/(sv@{j})*a@{j}(-1) + ((1-taup)*pi@{j}+(1-tau)*w*hc@{j}) - c@{j}; %overvej også -1 for de andre
@#endfor

    [name = 'Households Budget Constraint 24']
    a25 = (1+r)/sv24*a24(-1) + ((1-taup)*pi24+(1-tau)*w*hc24) + q25(+1) - c24; %overvej at fjerne +1 
    
[name = 'Bequest equation']
    q25 = (x74(-1)*fert_25(-74)*fert_25(-49)*n74(-1))/(n25); %fertilitet for bedsteforældre og forældre

@#for j in 25:J-32
    [name = 'Households Budget Constraint @{j}']
    a@{j+1} = (1+r)/sv@{j}*a@{j}(-1) + ((1-taup)*pi@{j}+(1-tau)*w*hc@{j}) - c@{j}; 
@#endfor

@#for j in 43:J-1
    [name = 'Households Budget Constraint @{j}']
    a@{j+1} = (1+r)/sv@{j}*a@{j}(-1) - c@{j} + d*w; 
@#endfor
    
    [name = 'Households Budget Constraint 74']
    a75 = (1+r)/sv74*a74(-1) - (fert_25(-73)*fert_25(-48)*x74) - c74 + d*w; %fertilitet for både bedsteforældre selv og deres børns fertilitet.

[name = 'Households Budget Constraint: terminal condition']
a75 = 0;

//Households: Financial Constraints (OBCs)
%--------------------------------------------------------------------------

@#for j in 1:J-32
    [name = 'Households Financial Constraint @{j}']
    min(lambda@{j},a@{j}+(D@{j}(+1)*opt_inc@{j}(+1))) = 0;
@#endfor
@#for j in 43:J
    [name = 'Households Financial Constraint @{j}']
    min(lambda@{j},a@{j}) = 0;
@#endfor

@#for j in 1:J-32
    opt_inc@{j} = w*hc@{j};
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
w  = (price*(((alpha)*(AK*K)^((sigma-1)/sigma)+(1-alpha)*(AL*L)^((sigma-1)/sigma))^(1/(sigma-1)))*(1-alpha)*AL^((sigma-1)/sigma)*L^(-1/sigma))/A_adj;
rk = (price*(((alpha)*(AK*K)^((sigma-1)/sigma)+(1-alpha)*(AL*L)^((sigma-1)/sigma))^(1/(sigma-1)))*(alpha)*AK^((sigma-1)/sigma)*K^(-1/sigma))/A_adj;             
Y  = (((alpha)*(AK*K)^((sigma-1)/sigma)+(1-alpha)*(AL*L)^((sigma-1)/sigma))^(sigma/(sigma-1)))/A_adj;
r = rk/e(-1) + (1-delta)*e/e(-1) - 1;
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

gov_debt*K=(gov_debt(-1)*K(-1)*(1+r(-1))+g(-1)*Y(-1)+ E(-1) -(gov_rev(-1)*(1-gov_deficit(-1)))); %ganger med K jf. dynare_ss_2024 og har tilføjet E(-1)
gov_deficit = (b(+1)*Y(+1)-gov_debt*K)/gov_rev;
gov_rev = (g*Y+E+r*gov_debt*K); %statens udgifter

tau*w*L = ((gov_rev)*(1-gov_deficit))-taup*PI;


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
    + n@{j}*c@{j}
@#endfor
);


[name = '(K) Kapital'] % Egg. bruger e(-1), men CG foreslår e. Egg. anvender a@{j}(-1), hvilket skyldes budgetbetingelsen.
K = (
@#for j in 1:J
   + n@{j}*a@{j}(-1)/sv@{j}
@#endfor
)/(e+gov_debt); %dividerer med gov_debt da K indgår i definitionen på gov_debt.


end;


%--------------------------------------------------------------------------
//TRANSITIONAL DYNAMICS BLOCK: Initial and Ending Values
%--------------------------------------------------------------------------

initval;

%----------------------------Endog.Variables-------------------------------

@#for j in 1:74
    n@{j}=ss_1970(@{j},1);
@#endfor
 
@#for j in 1:74
    c@{j}=c_scaled_1970(@{j},1);
@#endfor
 
@#for j in 1:75
    a@{j}=a_scaled_1970(@{j},1);
@#endfor
 
@#for j in 224:297
    lambda@{j-223}=ss_1970(@{j},1);
@#endfor
 
@#for j in 298:339
    pi@{j-297}=ss_1970(@{j},1);
@#endfor
 
x74=x74_scaled_1970(1,1); %disse er ændret
q25=ss_1970(341,1);
w=ss_1970(342,1);
rk=ss_1970(343,1);
r=ss_1970(344,1);
price=ss_1970(345,1);
tau=ss_1970(346,1);
PI=ss_1970(347,1);
Y=ss_1970(348,1);
N=ss_1970(349,1);
L=ss_1970(350,1);
C=ss_1970(351,1);
K=ss_1970(352,1);
gov_rev=ss_1970(353,1);
gov_deficit=ss_1970(354,1);
gov_debt=ss_1970(355,1);
E=ss_1970(361,1);
@#for j in 1:J-32
    opt_inc@{j} = ss_1970(342,1)*hc@{j};
@#endfor

%----------------------------Exog.Variables--------------------------------

fert_25=exo_matrix(1,1);
e=exo_matrix(1,2);
b=exo_matrix(1,3);
theta=exo_matrix(1,4);
AL= exo_matrix(1,5);
%g = exo_matrix(1,270); %hvis den skal være tidsvarierende


@#for j in 6:79
	s@{j-5}=exo_matrix(1,@{j});
@#endfor

@#for j in 80:153
	sv@{j-79}=exo_matrix(1,@{j});
@#endfor

@#for j in 154:227
	su@{j-153}=exo_matrix(1,@{j});
@#endfor

@#for j in 228:269
	D@{j-227}=exo_matrix(1,@{j});
@#endfor

end;

%options_.debug=1;
%resid;
%steady(maxit=1000, solve_algo=0);
%check;

%--------------------------------------------------------------------------

endval;

%----------------------------Endog.Variables-------------------------------

@#for j in 1:74
    n@{j}=ss_2024(@{j},1);
@#endfor

@#for j in 1:74
    c@{j}=c_scaled_2024(@{j},1);
@#endfor

@#for j in 1:75
    a@{j}=a_scaled_2024(@{j},1);
@#endfor

@#for j in 224:297
    lambda@{j-223}=ss_2024(@{j},1);
@#endfor
 
@#for j in 298:339
    pi@{j-297}=ss_2024(@{j},1);
@#endfor

x74=x74_scaled_2024(1,1); 
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
E = ss_2024(360);
@#for j in 1:J-32
    opt_inc@{j} = ss_2024(342,1)*hc@{j};
@#endfor

%----------------------------Exog.Variables--------------------------------

fert_25=exo_matrix(201,1); %erstatter 152 med 201
e=exo_matrix(201,2);
b=exo_matrix(201,3);
theta=exo_matrix(201,4);
AL=exo_matrix(201,5);
%g = exo_matrix(201,270); %hvis den skal være tidsvariende

%disse er ændret:
@#for j in 6:79
	s@{j-5}=exo_matrix(201,@{j});
@#endfor

@#for j in 80:153
	sv@{j-79}=exo_matrix(201,@{j});
@#endfor

@#for j in 154:227
	su@{j-153}=exo_matrix(201,@{j});
@#endfor

@#for j in 228:269
	D@{j-227}=exo_matrix(201,@{j});
@#endfor

end;
 
%--------------------------------------------------------------------------
//TRANSITIONAL DYNAMICS BLOCK: COMPUTATION
%--------------------------------------------------------------------------

perfect_foresight_setup(periods=199); %T-2
oo_.exo_simul=exo_matrix;
perfect_foresight_solver(stack_solve_algo=7, solve_algo=9);




%--------------------------------------------------------------------------
//SAVE FILES FOR COMPARISON WITH EGGERTSSON ET AL. (2019)
%--------------------------------------------------------------------------

save('dynare_transition.mat');

%--------------------------------------------------------------------------
//LATEX FILES
%--------------------------------------------------------------------------

%write_latex_dynamic_model;
%write_latex_static_model;
%write_latex_definitions;
%write_latex_parameter_table;
%collect_latex_files;
