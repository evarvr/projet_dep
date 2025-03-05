%clear;
close all;

global texp xexp tref xref E0 Es0 Ed0 tsw Itsw xsw kd modele k1F k1B k2F k2B k3F k3B k4 numfich

% Initialisation
modele = 'bbpp';  % Choix du modèle
E0 = 10;         % E0 fixé
Es0 = 0.61 * E0;
Ed0 = 0.39 * E0;
kd = 0.174;

% Lecture du fichier de données
numfich = input('Numero de fichier : ', 's');
filepath = ['data/mesures-' numfich '.dat'];
if exist(filepath, 'file') ~= 2
    error('Le fichier %s n''existe pas.', filepath);
end

% Chargement des données expérimentales
A = load(filepath);
texp = A(:,1);
xexp = A(:,2:end);
tsw = input(['  Temps de commutation du modele sur [0 ; ' num2str(texp(end)) ' ] : ']);

texp = texp(:); % S'assurer que texp est un vecteur colonne

% Paramètres d'optimisation
alpha = 10;
m = 140;
iteration_count = 0;
iterfailure = 0;

% Solution initiale

k1F = 5.6427; k1B = 4962.5543; k2F = 2103.0069; k2B = 456.3628;
k3F = 322.1195; k3B = 1214.6433; k4 = 1695.0416;

%k1F = 7.4912; k1B = 4963.6154; k2F = 2070.0574; k2B = 460.372; k3F = 320.1887; k3B = 1206.4745; k4 = 1716.7748;

%k1F = 0.0100; k1B = 5000.0000; k2F = 3333.3367; k2B = 1666.6733; k3F = 1666.6733; k3B = 5000.0000; k4 = 3333.3367;

% Meilleure solution trouvée
k1F_best = k1F; k1B_best = k1B; k2F_best = k2F; k2B_best = k2B;
k3F_best = k3F; k3B_best = k3B; k4_best = k4;

% Évaluation initiale
[error_stat, error_min, error_max, error_stat_rel, error_min_rel, error_max_rel] = ...
    simulate_and_evaluate(E0, texp, xexp, k1F, k1B, k2F, k2B, k3F, k3B, k4);

% Boucle d'optimisation
while iterfailure <= m
    iteration_count = iteration_count + 1;

    % Perturbation locale
    k1F = max(min(k1F + (rand() * 2 - 1) * alpha, 5000), 0);
    k1B = max(min(k1B + (rand() * 2 - 1) * alpha, 5000), 0);
    k2F = max(min(k2F + (rand() * 2 - 1) * alpha, 5000), 0);
    k2B = max(min(k2B + (rand() * 2 - 1) * alpha, 5000), 0);
    k3F = max(min(k3F + (rand() * 2 - 1) * alpha, 5000), 0);
    k3B = max(min(k3B + (rand() * 2 - 1) * alpha, 5000), 0);
    k4 = max(min(k4 + (rand() * 2 - 1) * alpha, 5000), 0);

    % Simulation avec les nouveaux paramètres
    [error_stat_new, error_min_new, error_max_new, error_stat_rel_new, error_min_rel_new, error_max_rel_new] = simulate_and_evaluate(E0, texp, xexp, k1F, k1B, k2F, k2B, k3F, k3B, k4);

    fprintf('\n=== ITERATION %d ===\n', iteration_count);
    fprintf('  Erreur stat = %.6f | Erreur min = %.6f | Erreur max = %.6f\n', ...
            error_stat_new, error_min_new, error_max_new);

    % Vérification de l'amélioration
    if error_stat_new < error_stat
        error_stat = error_stat_new;
        error_min = error_min_new;
        error_max = error_max_new;
        error_stat_rel = error_stat_rel_new;
        error_min_rel = error_min_rel_new;
        error_max_rel = error_max_rel_new;
        
        k1F_best = k1F;
        k1B_best = k1B;
        k2F_best = k2F;
        k2B_best = k2B;
        k3F_best = k3F;
        k3B_best = k3B;
        k4_best = k4;
        
        iterfailure = 0;
    else
        iterfailure = iterfailure + 1;
    end
end

% Affichage des résultats finaux
fprintf('\n=== MEILLEURE SOLUTION TROUVÉE ===\n');
fprintf('k1F = %.4f, k1B = %.4f, k2F = %.4f, k2B = %.4f, k3F = %.4f, k3B = %.4f, k4 = %.4f\n', ...
        k1F_best, k1B_best, k2F_best, k2B_best, k3F_best, k3B_best, k4_best);
fprintf('Erreur statistique = %.4f\n', error_stat);
fprintf('Erreur min = %.4f, Erreur max = %.4f\n', error_min, error_max);
fprintf('Erreur stat relative = %.4f\n', error_stat_rel);
fprintf('Erreur min relative = %.4f, Erreur max relative = %.4f\n', error_min_rel, error_max_rel);

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
