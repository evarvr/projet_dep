%clear all;
close all;

tic; % début du temps

global texp xexp tref xref E0 Es0 Ed0 tsw Itsw xsw kd modele k1F k1B k2F k2B k3F k3B k4 numfich

% Initialisation du programme
modele = 'bbpp'; % Choix du modèle
E0 = 10;         % E0 fixé
Es0 = 0.61 * E0; % Valeur dérivée d'E0
Ed0 = 0.39 * E0; % Valeur dérivée d'E0

% Lecture du numéro de fichier UNE SEULE FOIS
numfich = input('Numero de fichier : ', 's');
filepath = ['data/mesures-' numfich '.dat'];
if exist(filepath, 'file') ~= 2
    error('Le fichier %s n''existe pas.', filepath);
end

% Charger les données une seule fois
A = load(filepath);
texp = A(:, 1);
xexp = A(:, 2:end);
tsw = input(['  Temps de commutation du modele sur [0 ; ' num2str(texp(end)) ' ] : ']);

% Afficher les dimensions pour vérification
disp('Dimensions des données expérimentales :');
disp(size(xexp));

% Définir les plages des paramètres
param_ranges = [
    0.01, 5000;   % k1F
    0.01, 5000;   % k1B
    0.01, 5000;   % k2F
    0.01, 5000;   % k2B
    0.01, 5000;   % k3F
    0.01, 5000;   % k3B
    0.01, 5000    % k4
];

% Nombre d'échantillons
n_samples = 1000; % Nombre d'échantillons pour Latin Hypercube

% Générer des échantillons Latin Hypercube
param_space = lhsdesign(n_samples, 7); % 7 paramètres
param_values = zeros(n_samples, 7);
for i = 1:7
    param_values(:, i) = param_space(:, i) * (param_ranges(i, 2) - param_ranges(i, 1)) + param_ranges(i, 1);
end

% Initialisation des variables pour stocker les résultats
best_params_stat = NaN(1, 7); % Initialisation avec des valeurs NaN
best_params_min = NaN(1, 7);
best_params_max = NaN(1, 7);
best_error_stat = inf;  % Erreur initiale très élevée
best_error_max=inf;
best_error_min=inf;
best_error_min_rel=inf;
best_error_max_rel=inf;
best_error_stat_rel=inf;

% Boucle d'optimisation
disp('Démarrage de l''optimisation avec Latin Hypercube Sampling...');
for i = 1:n_samples
    % Extraire les paramètres pour l'itération actuelle
    k1F = param_values(i, 1);
    k1B = param_values(i, 2);
    k2F = param_values(i, 3);
    k2B = param_values(i, 4);
    k3F = param_values(i, 5);
    k3B = param_values(i, 6);
    k4 = param_values(i, 7);

    % Évaluation de l'erreur pour ces paramètres
    [error_stat, error_min, error_max, error_stat_rel, error_min_rel, error_max_rel] = simulate_and_evaluate(E0, texp, xexp, k1F, k1B, k2F, k2B, k3F, k3B, k4);

    % Mise à jour des meilleurs paramètres si une amélioration est trouvée
    if error_stat < best_error
        best_error_stat = error_stat;
        best_params_stat = [k1F, k1B, k2F, k2B, k3F, k3B, k4];
        best_error_stat_rel = error_stat_rel;  % Enregistrer l'erreur relative statistique
    end

    if error_min < best_error_min
        best_error_min = error_min;
        best_params_min = [k1F, k1B, k2F, k2B, k3F, k3B, k4];
        best_error_min_rel = error_min_rel;  % Enregistrer l'erreur relative minimale
    end

    if error_max < best_error_max
        best_error_max = error_max;
        best_params_max = [k1F, k1B, k2F, k2B, k3F, k3B, k4];
        best_error_max_rel = error_max_rel;  % Enregistrer l'erreur relative maximale
    end
end

% Affichage des résultats
%fprintf('Nombre total d''itérations : %d\n', iteration);

fprintf('Meilleurs paramètres (erreur statistique) : k1F = %.4f, k1B = %.4f, k2F = %.4f, k2B = %.4f, k3F = %.4f, k3B = %.4f, k4 = %.4f, Erreur = %.4f\n', ...
    best_params_stat, best_error_stat);

fprintf('Meilleurs paramètres (erreur minimale) : k1F = %.4f, k1B = %.4f, k2F = %.4f, k2B = %.4f, k3F = %.4f, k3B = %.4f, k4 = %.4f, Erreur = %.4f\n', ...
    best_params_min, best_error_min);

fprintf('Meilleurs paramètres (erreur maximale) : k1F = %.4f, k1B = %.4f, k2F = %.4f, k2B = %.4f, k3F = %.4f, k3B = %.4f, k4 = %.4f, Erreur = %.4f\n', ...
    best_params_max, best_error_max);

% Affichage des résultats
%fprintf('Nombre total d''itérations : %d\n', iteration);

fprintf('Meilleurs paramètres (erreur statistique) : k1F = %.4f, k1B = %.4f, k2F = %.4f, k2B = %.4f, k3F = %.4f, k3B = %.4f, k4 = %.4f, Erreur = %.4f\n', ...
    best_params_stat, best_error_stat);

fprintf('Meilleurs paramètres (erreur minimale) : k1F = %.4f, k1B = %.4f, k2F = %.4f, k2B = %.4f, k3F = %.4f, k3B = %.4f, k4 = %.4f, Erreur = %.4f\n', ...
    best_params_min, best_error_min);

fprintf('Meilleurs paramètres (erreur maximale) : k1F = %.4f, k1B = %.4f, k2F = %.4f, k2B = %.4f, k3F = %.4f, k3B = %.4f, k4 = %.4f, Erreur = %.4f\n', ...
    best_params_max, best_error_max);


% Affichage des erreurs relatives
fprintf('Erreur statistique relative : %.4f\n', best_error_stat_rel);
fprintf('Erreur minimale relative : %.4f\n', best_error_min_rel);
fprintf('Erreur maximale relative : %.4f\n', best_error_max_rel);



function [error_stat, error_min, error_max,error_stat_rel, error_min_rel ,error_max_rel] = simulate_and_evaluate(E0, texp, xexp, k1F, k1B, k2F, k2B, k3F, k3B, k4)
    % Appel du simulateur (Simulateur.m)
    run('trace_erreur.m');
    
        if exist('J1opt', 'var') && ~isempty(J1opt)
            error_stat = J1opt; % Utiliser l'erreur quadratique
            error_stat_rel=J_rel;
        else
            error_stat = inf;
            fprintf('J1opt n''est pas défini ou est vide. Erreur assignée à inf.\n');
        end   

        % Appeler une fonction objectif min-max
        [error_min,error_max, error_min_rel, error_max_rel] = objectif_minmax(E0, texp, xexp, k1F, k1B, k2F, k2B, k3F, k3B, k4);
   
end