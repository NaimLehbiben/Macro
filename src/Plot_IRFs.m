function [beta_T, beta_T_std] = Plot_IRFs(H, p, country_names, graph_title, Y, T, Y_control)

% --------------------------------------------------------------------------------------------------------------------------
% Plot_IRFs : Trace les fonctions de réponse impulsionnelle (IRFs) pour chaque pays.
%
% Entrées :
% - H : Horizon maximal auquel les IRFs sont calculées.
% - p : Nombre de retards pour les variables dépendantes dans la régression.
% - country_names : Tableau des noms des pays.
% - graph_title : Titre pour l'ensemble du graphique.
% - Y : Variable endogène.
% - T : Vecteur de l'anomalie de temperature.
% - Y_control : Matrice des variables de contrôle supplémentaires (optionnelle).
%
% Sorties :
% - beta_T : Matrice des coefficients estimés des IRFs pour chaque pays (H+1 lignes x nb_pays colonnes).
% - beta_T_std : Matrice des écarts-types des coefficients estimés des IRFs (H+1 lignes x nb_pays colonnes).
% --------------------------------------------------------------------------------------------------------------------


% Paramètres
nb_country = length(country_names);
beta_T = zeros(H+1,nb_country);
beta_T_std = zeros(H+1,nb_country);


% Boucle sur les pays
figure()
for i=1:nb_country

    % Calcul pour le mois m des prix au mois m-1
    Pm_m1 = lagmatrix(repmat(Y(:,i),1,H+1),1);

    % Calcul pour le mois m des prix aux horizons m,m+1, ...,m+H
    Pm_ph = lagmatrix(Y(:,i),0:-1:-H); 

    % Calcul des taux de croissance log
    DP = (log(Pm_ph) - log(Pm_m1))*100;

    % Elimination des données manquantes 
    idx_ = ~any(isnan(DP),2);
    DP = DP(idx_,:); % date deb : 02/1996 // date fin : 12/2021
    
    % Boucle sur les horizons et calcul des IFRs
    for h=1:H+1
    
        % Indice de départ
        ts = p+1;
    
        % Matrices de régression
        X = lagmatrix(DP(:,h),1:p); % lag de la variable 

         % Si Y_control est fourni, l'inclure dans la matrice X
        if ~isempty(Y_control)
            X = [T(ts:end,i) X(ts:end,:) Y_control(ts:end,:)];
        else
            X = [T(ts:end,i) X(ts:end,:)];
        end
    
        % Régression - HAC
        %[beta_cov, ~, reg_param] = hac(X, DP(ts:end,h), 'Display', 'off');
        %beta_T(h,i) = reg_param(1);
        %beta_T_std(h,i) = sqrt(beta_cov(1,1));

        % Régression - OLS classique (moins de temps)
        ols_mdl  = fitlm(X, DP(ts:end,h),'Intercept', false);

        % Beta (IFRs)
        reg_param = ols_mdl.Coefficients.Estimate;
        beta_T(h,i) = reg_param(1);

        % Std du beta (IFRs)
        beta_cov = NeweyWest_covmatrix(ols_mdl.Residuals.Raw,X);
        beta_T_std(h,i) = sqrt(beta_cov(1,1));
    
    end

    % Compteur
    disp(['Country ' num2str(i) ' ' country_names{i}])

    % Graphique
    subplot(3,4,i);
    x = (0:H)';

    % Intervalle de confiance à 68%
    beta_conf_inf = beta_T(:,i)-beta_T_std(:,i);
    beta_conf_sup = beta_T(:,i)+beta_T_std(:,i);
    fill([x; flipud(x)], [beta_conf_sup; flipud(beta_conf_inf)], ...
        'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');       

    % IRFs
    hold on
    line(xlim, [0 0], 'Color', 'k', 'LineWidth', 1);
    plot(0:H, beta_T(:,i), 'LineWidth',1.5);  
    hold off
    title([country_names{i}]);
    sgtitle(graph_title);
    
end


end