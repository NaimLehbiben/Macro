%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Projet Macroéconomie %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all;

%% Chargement des données

% Récupération des pays
opts = detectImportOptions("Data.xlsx");
preview("Data.xlsx", opts);
country_names = opts.VariableNames(2:end);
nb_country = length(country_names);

% Métrique climatique
anomaly = readmatrix("Data.xlsx","Sheet","Anomaly");
anomaly_precipitation = readmatrix("Data.xlsx","Sheet","Anomaly_precipitation");
anomaly_forecast = readmatrix("Data.xlsx","Sheet","Anomaly_2100_SSP2");
anomaly_forecast_SSP1 = readmatrix("Data.xlsx","Sheet","Anomaly_2100_SSP1");
anomaly_forecast_SSP3 = readmatrix("Data.xlsx","Sheet","Anomaly_2100_SSP3");
anomaly_forecast_SSP5 = readmatrix("Data.xlsx","Sheet","Anomaly_2100_SSP5");

% Metriques macroéconomique
hicp =  readmatrix("Data.xlsx","Sheet","HICP");
gdp =  readmatrix("Data.xlsx","Sheet","GDP");
dates = datetime(hicp(:,1), 'ConvertFrom', 'excel');


% Suppression dates
anomaly = anomaly(:,2:end);
anomaly_precipitation = anomaly_precipitation(:,2:end);
hicp = hicp(:,2:end);
gdp = gdp(:,2:end);

% SSP anomalie
anomaly_forecast = anomaly_forecast(2:end,2:end);
anomaly_forecast_SSP1 = anomaly_forecast_SSP1(2:end,2:end); % SSP1 2.6
anomaly_forecast_SSP3 = anomaly_forecast_SSP3(2:end,2:end);
anomaly_forecast_SSP5 = anomaly_forecast_SSP5(2:end,2:end);

%% 3.2.A - Map des anomalies

% Mapping Toolbox + borders package nécessaires

% % Colormap
% cmap = colormap(flipud(pink(100))); % specify the colormap
% % Données de temperétaure
% data = mean(T_m); % i'm assuming this is a placeholder
% datarange = [min(data) max(data)]; % nominal data range needs to be symmetric
% % Rescaling -> pour couleur
% coloridx = round((data-min(datarange))/range(datarange)*(size(cmap,1)-1) + 1);
% 
% % Map
% figure('Color','w')
% xh = worldmap('Europe');
% setm(xh,'MapProjection','miller')
% tightmap
% geoshow ('landareas.shp', 'FaceColor', 'white');
% for i = 1:1:nb_country
%     bordersm(country_names{i},'facecolor',cmap(coloridx(i),:)) 
% end
% colorbar
% colormap(cmap)
% clim(datarange) 
% 

%% 3.2.B - Map des anomalies forecastées
% 
% % Colormap
% cmap = colormap(flipud(pink(100))); % specify the colormap
% % Données de temperétaure
% data = mean(anomaly_forecast); % i'm assuming this is a placeholder
% datarange = [min(data) max(data)]; % nominal data range needs to be symmetric
% % Rescaling -> pour couleur
% coloridx = round((data-min(datarange))/range(datarange)*(size(cmap,1)-1) + 1);
% 
% % Map
% figure('Color','w')
% xh = worldmap('Europe');
% setm(xh,'MapProjection','miller')
% tightmap
% geoshow ('landareas.shp', 'FaceColor', 'white');
% for i = 1:1:nb_country
%     bordersm(country_names{i},'facecolor',cmap(coloridx(i),:)) 
% end
% colorbar
% colormap(cmap)
% clim(datarange) 

%% 3.3.A - Calcul des IRFs - HICP

% Paramètres
H_m = 24;
p_HICP = 6; % nombre de lag de l'HICP
L = size(anomaly,1) - p_HICP;
alpha = 0.90;

% Variables
T_m = anomaly(2:end,:); % % date deb : 02/1996 // date fin : 12/2021


%, Affichage des IRFs
[beta_T, beta_T_std] = Plot_IRFs(H_m, p_HICP, country_names, "IRFs - HICP", hicp, T_m, []);

% Ajustement de l'échelle des graphiques 
figHandles = findall(gcf, 'Type', 'axes');
for k = 1:length(figHandles)
    ylim(figHandles(k), [-0.2 0.2]);
end

%% 3.3.B - Calcul des IRFs - GDP

% Paramètres
H_q = 8;
p_gdp = 2; % nombre de lag du GDP
L = size(anomaly,1)/3 - p_gdp;
alpha = 0.90;

% Variables
T_q = zeros(size(anomaly,1)/3,nb_country);
for i=1:nb_country
    T_q(:,i) = sum(reshape(anomaly(:,i), 3, []));
end
T_q = T_q(2:end,:);


% Affichage des IRFs
[beta_gdp, beta_gqp_std] = Plot_IRFs(H_q, p_gdp, country_names, "IRFs - GDP", gdp, T_q, []);

% Ajustement de l'échelle des graphiques 
figHandles = findall(gcf, 'Type', 'axes');
for k = 1:length(figHandles)
    ylim(figHandles(k), [-0.5 0.7]);
    xlim(figHandles(k),[0 8]);
end
%% 4.1.A - Inflation projection 2100 - Equation 5 et Figure 8

% Betas
last_beta_T = beta_T(end,:);

% Beta du papier (Appendix table 1)
last_beta_T = [0.0069743, -0.0002493, -0.1510064, -0.0666084, ...
            -0.0636513, -0.0195247, -0.3789182, 0.0663101, ...
            -0.2787566, 0.0096166, -0.2413174, -0.2386734];

% Comparatif betas
disp('Betas - Temperature')
betas_names = {'Empirical betas', 'Paper betas'};
beta_T_table = table(beta_T(end,:)', last_beta_T', ...
    'VariableNames', betas_names, 'RowNames', country_names);
disp(beta_T_table)

% Prévision d'inflation 
forecast_inflation = anomaly_forecast.*last_beta_T; 
forecast_inflation_SSP1 = anomaly_forecast_SSP1.*last_beta_T; 
forecast_inflation_SSP5 = anomaly_forecast_SSP5.*last_beta_T; 

% Graphique (Figure 8)
figure()
for i=1:nb_country

    % Graphique par pays
    subplot(3,4,i);
    x=(2025:2100)';

    % IC - SSP1 et SSP5 (best and worst case scenarios)
    fill([x; flipud(x)], [forecast_inflation_SSP5(:,i); flipud(forecast_inflation_SSP1(:,i))], ...
        'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 

    % Trajectoires d'inflation SSP2
    hold on
    line(xlim, [0 0], 'Color', 'k', 'LineWidth', 1);
    plot(x, forecast_inflation(:,i), 'LineWidth',1.5);
    hold off

    title([country_names{i}]);
    ylim([-3 2]);
    xlim([2025 2100]);
    ax = gca; % Obtenir l'axe courant
    ax.XTick = min(x):25:max(x);

end

% Titre
sgtitle('Temperature-induced changes in inflation (HICP) for each euro area economy (2025-2100)');

%% 4.1.B - GDP projection 2100 - Equation 6 et Figure 9

% Betas
last_beta_gdp = beta_gdp(end,:);

% Beta du papier (Appendix table 1)
last_beta_gdp =[-0.4695802, -0.2587495, -0.4326935, -0.2465675, ...
           -0.2011965, -0.3843606, -0.4510485, -0.1377134, ...
           -0.3133269, -0.1869576, 0.1762825, 0.1618944] ;

% Comparatif betas
disp('Betas - GDP')
betas_names = {'Empirical betas', 'Paper betas'};
beta_gdp_table = table(beta_gdp(end,:)', last_beta_gdp', ...
    'VariableNames', betas_names, 'RowNames', country_names);
disp(beta_gdp_table);

% Prévision de GDP
forecast_gdp = anomaly_forecast.* last_beta_gdp; 
forecast_gdp_SSP1 = anomaly_forecast_SSP1.* last_beta_gdp; 
forecast_gdp_SSP5 = anomaly_forecast_SSP5.* last_beta_gdp; 

% Graphique
figure()
for i=1:nb_country

    % Graphique du pays
    subplot(3,4,i);
    x=(2025:2100)';

    % IC - SSP1 et SSP5 (best and worst case scenarios)
    fill([x; flipud(x)], [forecast_gdp_SSP5(:,i); flipud(forecast_gdp_SSP1(:,i))], ...
        'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Trajectoires PIB SSP2
    hold on
    line(xlim, [0 0], 'Color', 'k', 'LineWidth', 1);
    plot(2025:2100, forecast_gdp(:,i), 'LineWidth',1.5);
    hold off

    title([country_names{i}]);
    ylim([-3 2]);
    xlim([2025 2100]);
    ax = gca; % Obtenir l'axe courant
    ax.XTick = min(x):25:max(x);

end

% Titre
sgtitle('Temperature-induced changes in GDP per capita for each euro area economy (2025-2100)');

%% 4.2 - Temperature induced policy rate - Figure 10

% Equation 8
nb_years_forecast = length(anomaly_forecast);
rate_forecast = ones(nb_years_forecast,nb_country)+1.5*forecast_inflation+0.5*forecast_gdp;
rate_forecast_SSP1 = ones(nb_years_forecast,nb_country)+1.5*forecast_inflation_SSP1+0.5*forecast_gdp_SSP1;
rate_forecast_SSP5 = ones(nb_years_forecast,nb_country)+1.5*forecast_inflation_SSP5+0.5*forecast_gdp_SSP5;

% Graphique
figure()
for i=1:nb_country

    % Graphique par pays
    subplot(3,4,i);
    x = (2025:2100)';

    % IC - SSP1 et SSP5 (best and worst case scenarios)
    fill([x; flipud(x)], [rate_forecast_SSP5(:,i); flipud(rate_forecast_SSP1(:,i))], ...
        'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');  

    % Graphique des politiques monétaires induite par les anomalies
    hold on
    line(xlim, [1 1], 'Color', 'k', 'LineWidth', 1);
    plot(x, rate_forecast(:,i), 'LineWidth',1.5);
    hold off

    title([country_names{i}]);
    xlim([2025 2100]);
    ylim([-5 2]);
    xlim([2025 2100]);
    ax = gca; % Obtenir l'axe courant
    ax.XTick = min(x):25:max(x);

end

% Titre
sgtitle('Temperature-induced changes in policy rate for each euro area economy (2025-2100)');

%% 4.3.A - Monetary policy stress - Equation 9 et Figure 5

% Matrice monetary stress
monetary_stress = zeros(nb_years_forecast,nb_country);
monetary_stress_SSP1 = zeros(nb_years_forecast,nb_country);
monetary_stress_SSP5 = zeros(nb_years_forecast,nb_country);

% Boucle sur les pays
figure()
for i=1:nb_country

    % Monetary stress SSP2
    rate_other_countries = rate_forecast(:, [1:i-1, i+1:end]); % Matrice des forecast de taux d'intérêt en excluant colonne du pays i
    mean_rate_other_countries = mean(rate_other_countries,2); % Moyenne des taux forecast sans ce pays à chaque date
    monetary_stress(:,i)=rate_forecast(:,i)-mean_rate_other_countries; % Monetary stress pour ce pays

    % Monetary stress SSP1
    rate_other_countries = rate_forecast_SSP1(:, [1:i-1, i+1:end]); % Matrice des forecast de taux d'intérêt en excluant colonne du pays i
    mean_rate_other_countries = mean(rate_other_countries,2); % Moyenne des taux forecast sans ce pays à chaque date
    monetary_stress_SSP1(:,i)=rate_forecast_SSP1(:,i)-mean_rate_other_countries; % Monetary stress pour ce pays

    % Monetary stress SSP5
    rate_other_countries = rate_forecast_SSP5(:, [1:i-1, i+1:end]); % Matrice des forecast de taux d'intérêt en excluant colonne du pays i
    mean_rate_other_countries = mean(rate_other_countries,2); % Moyenne des taux forecast sans ce pays à chaque date
    monetary_stress_SSP5(:,i)=rate_forecast_SSP5(:,i)-mean_rate_other_countries; % Monetary stress pour ce pays

    % Affichage graphique
    subplot(3,4,i);
    x = (2025:2100)';
   
    % IC - SSP1 et SSP5
    fill([x; flipud(x)], [monetary_stress_SSP5(:,i); flipud(monetary_stress_SSP1(:,i))], ...
        'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 

    hold on
    title([country_names{i}]);
    line(xlim, [0 0], 'Color', 'k', 'LineWidth', 1);
    plot(x, monetary_stress(:,i), 'LineWidth',1.5);
    hold off
    
    xlim([2025 2100]);
    ylim([-4 2]);
    xlim([2025 2100]);
    ax = gca; % Obtenir l'axe courant
    ax.XTick = min(x):25:max(x);
end

% Titre
sgtitle('Monetary stress of euro area economies (2025-2100)');

%% 4.3.B - GDP weighted aggregate monetary stress - Figure 7

% Poids GDP par rapport au GDP de toute l'euro area 2021
gdp_weights = readmatrix("Data.xlsx","Sheet","GDP_weights");
gdp_weights_avg = gdp_weights(2:end);

% SSP2
monatery_stress_abs = abs(monetary_stress);% Valeurs absolues monetary stress
monatery_stress_aggregate = monatery_stress_abs*gdp_weights_avg';% GDP-weighted aggregate monetary stress 

% SSP1
monatery_stress_abs_SSP1 = abs(monetary_stress_SSP1);% Valeurs absolues monetary stress
monatery_stress_aggregate_SSP1 = monatery_stress_abs_SSP1*gdp_weights_avg';% GDP-weighted aggregate monetary stress 

% SSP5
monatery_stress_abs_SSP5 = abs(monetary_stress_SSP5);% Valeurs absolues monetary stress
monatery_stress_aggregate_SSP5 = monatery_stress_abs_SSP5*gdp_weights_avg';% GDP-weighted aggregate monetary stress 

% Graphique
figure()
x = (2025:2100)';

% IC - SSP1/SSP5
fill([x; flipud(x)], [monatery_stress_aggregate_SSP5; flipud(monatery_stress_aggregate_SSP1)], ...
    'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on;
line(xlim, [0 0], 'Color', 'k', 'LineWidth', 1);
plot(x, monatery_stress_aggregate, 'LineWidth',1.5);
title("GDP-weighted aggregate monetary stress (2025-2100)");
xlim([2025 2100]);
ax = gca; % Obtenir l'axe courant
ax.XTick = min(x):25:max(x);
hold off;

%% 5.1.A - Average temperature increase 2025-2100

% Les chiffres donnés dans 5.1 : anomaly_2100 - anomaly_2025

disp("Augmentation moyenne de la température entre 2025 et 2100 selon le scénario : ");
mean_anomaly_SSP1 = mean(anomaly_forecast_SSP1(end,:)-anomaly_forecast_SSP1(1,:));
mean_anomaly_SSP2 = mean(anomaly_forecast(end,:)-anomaly_forecast(1,:));
mean_anomaly_SSP3 = mean(anomaly_forecast_SSP3(end,:)-anomaly_forecast_SSP3(1,:));
mean_anomaly_SSP5 = mean(anomaly_forecast_SSP5(end,:)-anomaly_forecast_SSP5(1,:));
fprintf('SSP1 : %3f\n', mean_anomaly_SSP1);
fprintf('SSP2 : %3f\n', mean_anomaly_SSP2);
fprintf('SSP3 : %3f\n', mean_anomaly_SSP3);
fprintf('SSP5 : %3f\n', mean_anomaly_SSP5);

%% 5.1.B - Projected changes in temperatures Figure 11

% Anomalie 2100 - 2025, pour chaque scénario
data = [anomaly_forecast_SSP1(end,:)-anomaly_forecast_SSP1(1,:); 
        anomaly_forecast(end,:)-anomaly_forecast(1,:);
        anomaly_forecast_SSP3(end,:)-anomaly_forecast_SSP3(1,:);
        anomaly_forecast_SSP5(end,:)-anomaly_forecast_SSP5(1,:)]'; 

% Graphique histogramme
figure;
h = bar(data, 'grouped');

% Couleurs des barres --> dégradé de bleu
colors = [
    0.728, 0.867, 0.902; % Bleu clair
    0.429, 0.808, 0.980;  % Bleu ciel
    0, 0, 1;             % Bleu
    0, 0, 0.545         % Bleu foncé 
    ];

for k = 1:length(h)
    h(k).FaceColor = colors(k, :);
end

% Titres graphique et axes
ylabel('Temperature Anomaly');
title("Projected changes in temperatures (2025-2100)");
set(gca, 'XTickLabel', country_names); % noms pays axe abscisses

% Légende
legend({'SSP1', 'SSP2', 'SSP3', 'SSP5'}, 'Location', 'northeastoutside');

% Ajuster la mise en page pour que tout soit lisible
set(gca, 'FontSize', 12);
grid on;


%% Calcul des IRFs - HICP - Précipitations

% Paramètres
H_m = 24;
p_HICP = 6; % nombre de lag de l'HICP
L = size(anomaly_precipitation,1) - p_HICP;
alpha = 0.90;

% Variables
T_m = anomaly_precipitation(2:end,:); % % date deb : 02/1996 // date fin : 12/2021

graph_title = 'Figure 2: Effect of temperature anomaly precipitation on inflation rate (1996m1, 2021m12)';

% Affichage des IRFs
[beta_T_p, beta_T_p_std] = Plot_IRFs(H_m, p_HICP, country_names, graph_title, hicp, T_m, []);


% Ajustement de l'échelle des graphiques 
figHandles = findall(gcf, 'Type', 'axes');
for k = 1:length(figHandles)
    ylim(figHandles(k), [-0.008 0.008]);
end

%% Calcul des IRFs - GDP - Température -- Précipitations

% Paramètres
H_q = 8;
p_GDP = 2; 
L = size(anomaly_precipitation,1) - p_GDP;
alpha = 0.90;

% Variables
T_q_p = zeros(size(anomaly_precipitation,1)/3,nb_country);
for i=1:nb_country
    T_q_p(:,i) = sum(reshape(anomaly_precipitation(:,i), 3, []));
end
T_q_p = T_q_p(2:end,:);

graph_title = 'Figure 3:  Effect of temperature anomaly precipitation on GDP per capita (1996Q1, 2021Q4)';

% Affichage des IRFs
[beta_gpd_p, beta_gpd_p_std] = Plot_IRFs(H_q, p_GDP, country_names, graph_title, gdp, T_q_p, []);

% Ajustement de l'échelle des graphiques 
figHandles = findall(gcf, 'Type', 'axes');
for k = 1:length(figHandles)
    ylim(figHandles(k), [-0.05 0.1]);
    xlim(figHandles(k),[0 8]);
end

%% Calcul des IRFs - HICP - With Control Variables

% Paramètres
H_m = 24;
p_HICP = 6; % nombre de lag de l'HICP
L = size(anomaly,1) - p_HICP;
alpha = 0.90;

% Variables
T_m = anomaly(2:end,:); % % date deb : 02/1996 // date fin : 12/2021


% --- Variables de controle  --- %

% Indicatrice saisonnière
months = repmat(1:12,1,size(anomaly,1)/12)';
seasonal_dummies = dummyvar(categorical(months));
%seasonal_dummies = seasonal_dummies(:,1:end-1);

% Trend 
trend = (1 : 1/12 : (1 + (size(anomaly,1)-1)*1/12))';
trend_quad = trend.^2;
trend_cub = trend.^3;
Y_control = [seasonal_dummies(2:end, :), trend(2:end, :), trend_quad(2:end, :), trend_cub(2:end, :)];

%, Affichage des IRFs
[beta_T, beta_T_std] = Plot_IRFs(H_m, p_HICP, country_names, "IRFs - HICP", hicp, T_m, Y_control);

% Ajustement de l'échelle des graphiques 
figHandles = findall(gcf, 'Type', 'axes');
for k = 1:length(figHandles)
    ylim(figHandles(k), [-0.2 0.2]);
end

%% Calcul des IRFs - GDP - With Control Variables

% Paramètres
H_q = 8;
p_gdp = 2; % nombre de lag du GDP
L = size(anomaly,1)/3 - p_gdp;
alpha = 0.90;

% Variables
T_q = zeros(size(anomaly,1)/3,nb_country);
for i=1:nb_country
    T_q(:,i) = sum(reshape(anomaly(:,i), 3, []));
end
T_q = T_q(2:end,:);


% --- Variables de controle  --- %

% Indicatrice saisonnière
quarters = repmat(1:4,1,ceil(size(anomaly,1)/12))';
seasonal_dummies = dummyvar(categorical(quarters));

% Trend 
trend = (1 : 1/4 : (1 + (size(T_q,1))*1/4))';
trend_quad = trend.^2;
trend_cub = trend.^3;
Y_control = [seasonal_dummies(2:end, :) trend(2:end, :), trend_quad(2:end, :), trend_cub(2:end, :)];

% Affichage des IRFs
[beta_gdp, beta_gqp_std] = Plot_IRFs(H_q, p_gdp, country_names, "IRFs - GDP", gdp, T_q, Y_control);

% Ajustement de l'échelle des graphiques 
figHandles = findall(gcf, 'Type', 'axes');
for k = 1:length(figHandles)
    ylim(figHandles(k), [-0.7, 0.7]);
    xlim(figHandles(k),[0 8]);
end

