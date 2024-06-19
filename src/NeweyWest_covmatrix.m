function V_nw = NeweyWest_covmatrix(e,X,H)

% Fonction calculant la matrice de covariance robuste à l'autocorrélation
% et l'hétéroscédasticité de Newes-West
% -------------------------------------------------------------------------
% Inputs:
%--------------------------------------------------------------------------
% - e : vecteur des résidus (N,1)
% - X : matrice des variables de la régression (N,k)
% - H : retard max de l'autocorrélation des résidus
%--------------------------------------------------------------------------
% Output:
%--------------------------------------------------------------------------
% - V_nw = matrice de variance covariance robuste de Newey-West

%% Paramètres

% Taille de l'échantillon
[N,k] = size(X);

% Si aucun retard spécifié: automatic bandwidth selection (Newey & West)
if nargin < 3
    H = floor(12*((N/100)^(2/9))); 
    % Autres coefficients possibles : 12, ceil(4*((N/100)^(2/9)));
end    


%% Calibration de la matrice de covariance robuste

% Matrice des erreurs x variables
Xe = X.*e;

% Sigma 0 
Q = Xe(1:N,:)'*Xe(1:N,:);

% Boucle sur les retards h
for h = 1:H
    % Coefficient de Bartlett
    w_h = 1 - h/(H+1);
    % Echantillon d'étude
    t = h+1:N;
    % Incrément de la matrice de covariance avec omega h
    Q = Q + w_h*(Xe(t,:)'*Xe(t-h,:) + Xe(t-h,:)'*Xe(t,:));
end

% Matrice de covariance robuste des résidus
Q = Q./(N-k);

% Matrice de covariance robuste des paramètres
V_nw = N.*((X'*X)\Q/(X'*X));

end