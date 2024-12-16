 
%%Individual Variables [ps_full in Eggertsson et al. (2019)]%%%%%%%%%%%%%%%
T=201; %ændr indeks nedenunder
%Population of each generation (n{j})
pop_indiv_dyn=Simulated_time_series.data(1:T,1:74);
pop_indiv_egg=ps_full.demog.pop;
%Consumption of each generation (c{j})
c_indiv_dyn=Simulated_time_series.data(1:T,75:148);
c_indiv_egg=ps_full.opt.C;
%Asset of each generation (a{j})
a_indiv_dyn=Simulated_time_series.data(1:T,149:222);
a_indiv_egg=ps_full.opt.a;
%Profit of each generation (p{j})
profit_indiv_dyn=Simulated_time_series.data(1:T,298:339);
profit_indiv_egg=ps_full.opt.profit;
%Bequest received by generation 32 (q32)
br_indiv_dyn=Simulated_time_series.data(1:T,341);
br_indiv_egg=ps_full.opt.br;
%Bequest given by generation 56 (x56)
bgo_indiv_dyn=Simulated_time_series.data(1:T,340);
bgo_indiv_egg=ps_full.opt.bgo;

%%Economy Variables [economy_full in Eggertsson et al. (2019)]%%%%%%%%%%%%%

%Total Kapital (K)
K_dyn=Simulated_time_series.data(1:T,352);
K_egg=economy_full.ag.K;
%Total Labor (L)
L_dyn=Simulated_time_series.data(1:T,350);
L_egg=economy_full.ag.L;
%Total Income (Y)
Y_dyn=Simulated_time_series.data(1:T,348);
Y_egg=economy_full.ag.Y;
%Total Consumption (C)
C_dyn=Simulated_time_series.data(1:T,351);
C_egg=economy_full.ag.C;
%Total Profit (PI)
PI_dyn=Simulated_time_series.data(1:T,347);
PI_egg=economy_full.ag.profit;
%Total Population (N)
N_dyn=Simulated_time_series.data(1:T,349);
N_egg=economy_full.ag.pop;
%Y/N og Y/L
y_capita = Y_dyn ./ N_dyn;  % Element-wise division
y_worker = Y_dyn ./ L_dyn;  % Element-wise division

%%Government Variables [gov_full in Eggertsson et al. (2019)]%%%%%%%%%%%%%%

%Government Debt
gov_debt_dyn=Simulated_time_series.data(1:T,355);
gov_debt_egg=gov_full.debt;
%Government Deficit
gov_deficit_dyn=Simulated_time_series.data(1:T,354);
gov_deficit_egg=gov_full.deficit;
%Government Revenues 
gov_tax_revt_dyn=Simulated_time_series.data(1:T,353);
gov_tax_revt_egg=gov_full.tax.revt;
%Government Taxation (tau)
gov_tax_dyn=Simulated_time_series.data(1:T,346);
gov_tax_egg=ps_full.tax.wa(:,1);

%%Prices Variables [prices_full in Eggertsson et al. (2019)]%%%%%%%%%%%%%%%

%Rental k (rk)
rentk_dyn=Simulated_time_series.data(1:T,343);
rentk_egg=prices_full.rentk;
%Interest Rate (r)
r_dyn=Simulated_time_series.data(1:T,344);
r_egg=prices_full.r;
%Wages (w)
wages_dyn=Simulated_time_series.data(1:T,342);
wages_egg=prices_full.wages;

%%PLOT COMPARISON OF MAIN VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Individual Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaultfigurecolor',[1 1 1 ])
set( gca                       , ...
    'FontName'   , 'Arial' );

%Population of each generation (n(j))
figure;
for j=1:56
    %figure;
    plot(pop_indiv_egg(:,j))
    hold on
    plot(pop_indiv_dyn(:,j))
    hold on
end
title('n(j)')

%Consumption of each generation (c(j))
figure;
for j=1:56
    %figure;
    plot(c_indiv_egg(:,j))
    hold on
    plot(c_indiv_dyn(:,j))
    hold on
end
title('c(j)')

%Asset of each generation (a(j))
figure;
for j=1:56
    %figure;
    plot(a_indiv_egg(2:152,j))
    hold on
    plot(a_indiv_dyn(1:151,j))
    hold on
end
title('a(j)')

%Profit of each generation (pi(j))
figure;
for j=1:40
    plot(profit_indiv_egg(:,j))
    hold on
    plot(profit_indiv_dyn(:,j))
    hold on
end
title('pi(j)')

%q32 and ps_full.opt.br(1,32)
figure;
plot(br_indiv_egg(1:150,32),'--b','LineWidth', 2.5);
hold on
plot(br_indiv_dyn(1:200),'-r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title({'Bequest received, q'})
%xlim([t(1) t(end)])
xtickangle(45);
set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );


saveas(gcf,'Figures/Bequest','epsc')


%x56 and ps_full.opt.bgo(1,56)
figure;
plot(bgo_indiv_egg(:,56),'--b','LineWidth', 2.5);
hold on
plot(bgo_indiv_dyn,'-r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title({'Bequest given, x'})
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

saveas(gcf,'Figures/x56','epsc')


%%%%Economy Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Aggregate Capital (K)
figure;
plot(K_egg,'--b','LineWidth', 2.5);
hold on
plot(K_dyn,'-r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Aggregate Capital (K)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );


saveas(gcf,'Figures/Aggregate_Kapital','epsc')


%Aggregate Labor (L)
figure;
plot(L_egg,'--b','LineWidth', 2.5);
hold on
plot(L_dyn,'-r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Aggregate Labor (L)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );


saveas(gcf,'Figures/Aggregate_Labor','epsc')


%Aggregate Income (Y)
figure;
plot(Y_egg,'--b','LineWidth', 2.5);
hold on
plot(Y_dyn,'-r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Aggregate Income (Y)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

saveas(gcf,'Figures/Aggregate_Income','epsc')

%Aggregate Consumption (C)
figure;
plot(C_egg,'--b','LineWidth', 2.5);
hold on
plot(C_dyn,'-r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Aggregate Consumption (C)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

saveas(gcf,'Figures/Aggregate_Consumption','epsc')


%Aggregate Profit (PI)
figure;
plot(PI_egg,'--b','LineWidth', 2.5);
hold on
plot(PI_dyn,'-r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Aggregate Profit (PI)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );


saveas(gcf,'Figures/Aggregate_Profit','epsc')


%Aggregate Population (N)
figure;
plot(N_egg,'--b','LineWidth', 2.5);
hold on
plot(N_dyn,'-r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Aggregate Population (N)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

saveas(gcf,'Figures/Aggregate_Population','epsc')

%BNP per capita (Y/N)
figure;
plot(y_capita, '-r', 'LineWidth', 2.5); hold on;  % Solid red line for y_capita
plot(y_worker, '--b', 'LineWidth', 2.5);          % Solid blue line for y_worker
legend('Y/N, Y/L')
title('BNP per capita og arbejder')

%xlim([t(1) t(end)]) % Uncomment this and adjust limits if you want to set custom limits

xtickangle(45);  % Rotating the x-axis labels

% Customizing the appearance of the axes
set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'     , ...  % No grid lines
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0       );


% Save the figure
saveas(gcf, 'Figures/BNP_per_capita_og_arbejder', 'epsc')



%%Government Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Government Debt
figure;
plot(gov_debt_egg,'--b','LineWidth', 2.5);
hold on
plot(gov_debt_dyn,'-r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Government Debt')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

saveas(gcf,'Figures/Gov_debt','epsc')


%Government Deficit
figure;
plot(gov_deficit_egg,'--b','LineWidth', 2.5);
hold on
plot(gov_deficit_dyn,'-r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Government Deficit')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

saveas(gcf,'Figures/Gov_deficit','epsc')


%Government Revenues
figure;
plot(gov_tax_revt_egg,'--b','LineWidth', 2.5);
hold on
plot(gov_tax_revt_dyn,'-r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Government Tax Revenues')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

saveas(gcf,'Figures/Gov_revenues','epsc')


%Government Taxation (tau)
figure;
plot(gov_tax_egg,'--b','LineWidth', 2.5);
hold on
plot(gov_tax_dyn,'-r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Government Taxation (tau)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

saveas(gcf,'Figures/Gov_Tax','epsc')


%%%%Prices Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Rental k (rentk)
figure;
plot(rentk_egg,'--b','LineWidth', 2.5);
hold on
plot(rentk_dyn,'-r','LineWidth', 2.5);
legend('Eggertsson et al. (2019)', 'Dynare')
title('Rental k (rk)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );


saveas(gcf,'Figures/rk','epsc')


%Interest rate (r)
figure;
plot(r_egg,'--b','LineWidth', 2.5);
hold on
plot(r_dyn,'-r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Interest Rate (r)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

saveas(gcf,'Figures/r','epsc')

%Wage (w)
figure;
plot(wages_egg,'--b','LineWidth', 2.5);
hold on
plot(wages_dyn,'-r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Wage (w)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

saveas(gcf,'Figures/wage','epsc')

%% Figure 8. Transition Path of the Natural Rate of Interest (pag. 42 paper)

r_egg_cut=r_egg(1:61);
r_dyn_cut=r_dyn(1:61);

[Trend_r_egg,Cyclical_r_egg] = hpfilter(r_egg_cut, 6.25);
%[Trend_r_egg,Cyclical_r_egg] = one_sided_hp_filter_kalman(r_egg_cut, 6.25);
[Trend_r_dyn,Cyclical_r_dyn] = hpfilter(r_dyn_cut, 6.25);
%[Trend_r_dyn,Cyclical_r_dyn] = one_sided_hp_filter_kalman(r_dyn_cut, 6.25);

%HP Trend Interest rate (r)
year_vec = [1:61]' + 1969;


figure;
plot(year_vec, Trend_r_egg,'--b','LineWidth',2);
%plot(year_vec, r_egg_cut,'LineWidth',2);

hold on
plot(year_vec, Trend_r_dyn,'-r','LineWidth',2);
%plot(year_vec, r_dyn_cut,'LineWidth',2);

hold on;
yline(0); hold on;
%legend('Eggertsson et al. (2019)', 'Dynare');
xlabel('Years');
%title('HP Trend Interest Rate % (r)');
legend('Eggertsson et al. (2019)', 'Dynare')

xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

saveas(gcf,'Figures/int_rate_hp','epsc')
