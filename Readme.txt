
Dansk:
Alle filerne og al data påkrevet for at simulere modellen (baseline) ligger i denne mappe.
Modellen simuleres ved at åbne "Main.m" og køre det.
Dette script aktiverer også dynare_ss_1970, dynare_ss_2024 og dynare_transition.
Hvis modellen skal rekalibreres, anvendes dynare_ss_1970_calibr, dynare_ss_2024_calibr, calibration_ss_1970 og calibration_ss_2024 også.

1970 henviser til startåret, som er 1970, hvor modellen starter i steady state. 
Modellen slutter i steady state 200 år senere, i 2170. Men da de eksogene variable holdes faste efter 2024, kaldes slutåret får 2024 i navngivningen af filerne. 

Undermapperne der starter med "+dynare" og "dynare" anvendes ikke.
De anvendte figurer i specialet ligger i mappen "plots". 
Nogle af dem vil ikke være simuleret i baseline, men kommer fra alternative simuleringer.

Data for de eksogene variable til at simulere modellen ligger i undermappen "data". 
De oprindelige data samt andre data ligger i mappen "Ny data" udenfor.
Alle eksogene variable (bort set fra hc), der anvendes i transitionen, ligger i "exo_matrix_start100".

English:
All the files and data required to simulate the model (baseline) are located in this folder.  
The model can be simulated by opening "Main.m" and running it.  
This script also executes `dynare_ss_1970`, `dynare_ss_2024`, and `dynare_transition`.  
If the model needs to be recalibrated, `dynare_ss_1970_calibr`, `dynare_ss_2024_calibr`, `calibration_ss_1970`, and `calibration_ss_2024` are also used.  

The year 1970 refers to the starting year, where the model begins in steady state.  
The model ends in steady state 200 years later, in 2170. However, since the exogenous variables remain constant after 2024, the terminal year is referred to as 2024 in the naming of the files.  

The subfolders starting with "+dynare" and "dynare" are not used.  
The figures used in the thesis are located in the "plots" folder. Some of these figures are not simulated in the baseline but come from alternative simulations.  

Data for the exogenous variables used to simulate the model is located in the "data" subfolder.  
The original data, along with other datasets, can be found in the folder "Ny data" outside this folder. 
Alle exogenous variables (except hc) used in the transition lie in "exo_matrix_start100".
 