clear;
close all;

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
alpha = 500;
m = 4;
iteration_count = 0;
iterfailure = 0;

% Solution initiale
%X = [11.7438, 3201.0492, 4015.9129, 1225.6538, 320.6112, 1315.7416, 513.6018];
X=[22.1147, 4671.6712, 3657.6492, 4105.9314, 706.3910, 3507.4235, 3490.8227];

% Extraction des paramètres
k1F = X(1); k1B = X(2); k2F = X(3); k2B = X(4);
k3F = X(5); k3B = X(6); k4 = X(7);

% Évaluation initiale
[error_stat, error_min, error_max, error_stat_rel, error_min_rel, error_max_rel] = ...
    simulate_and_evaluate(E0, texp, xexp, tsw, k1F, k1B, k2F, k2B, k3F, k3B, k4);

% Boucle d'optimisation
while iterfailure <= m
    iteration_count = iteration_count + 1;

    % Génération d'une nouvelle solution
    Xp = X + (rand(1, 7) * 2 - 1) * alpha;
    Xp = max(min(Xp, 5000), 0);

    % Extraction des nouveaux paramètres
    k1Fp = Xp(1); k1Bp = Xp(2); k2Fp = Xp(3); k2Bp = Xp(4);
    k3Fp = Xp(5); k3Bp = Xp(6); k4p = Xp(7);

    % Vérification des valeurs envoyées à trace_erreur2
    fprintf('\n=== ITERATION %d ===\n', iteration_count);
    fprintf('Envoi à trace_erreur2: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', ...
            k1Fp, k1Bp, k2Fp, k2Bp, k3Fp, k3Bp, k4p);

    % Simulation avec les nouveaux paramètres
    [error_stat_Xp, error_min_Xp, error_max_Xp, error_stat_rel_Xp, error_min_rel_Xp, error_max_rel_Xp] = ...
        simulate_and_evaluate(E0, texp, xexp, tsw, k1Fp, k1Bp, k2Fp, k2Bp, k3Fp, k3Bp, k4p);

    fprintf('  Erreur stat = %.6f | Erreur min = %.6f | Erreur max = %.6f\n', ...
            error_stat_Xp, error_min_Xp, error_max_Xp);

    % Vérification de l'amélioration
    if error_stat_Xp < error_stat
        X = Xp;
        k1F = X(1); k1B = X(2); k2F = X(3); k2B = X(4);
        k3F = X(5); k3B = X(6); k4 = X(7);

        error_stat = error_stat_Xp;
        error_min = error_min_Xp;
        error_max = error_max_Xp;
        error_stat_rel = error_stat_rel_Xp;
        error_min_rel = error_min_rel_Xp;
        error_max_rel = error_max_rel_Xp;

        iterfailure = 0;
    else
        iterfailure = iterfailure + 1;
    end

    % Pause pour analyse si nécessaire (décommente pour activer)
    % pause;
end

% Affichage des résultats finaux
fprintf('\n=== MEILLEURE SOLUTION TROUVÉE ===\n');
fprintf('k1F = %.4f, k1B = %.4f, k2F = %.4f, k2B = %.4f, k3F = %.4f, k3B = %.4f, k4 = %.4f\n', ...
        X(1), X(2), X(3), X(4), X(5), X(6), X(7));

fprintf('Erreur statistique = %.4f\n', error_stat);
fprintf('Erreur min = %.4f, Erreur max = %.4f\n', error_min, error_max);
fprintf('Erreur stat relative = %.4f\n', error_stat_rel);
fprintf('Erreur min relative = %.4f, Erreur max relative = %.4f\n', error_min_rel, error_max_rel);

%% **Fonction mise à jour**
function [error_stat, error_min, error_max, error_stat_rel, error_min_rel, error_max_rel] = simulate_and_evaluate(E0, texp, xexp, tsw, k1F, k1B, k2F, k2B, k3F, k3B, k4)

    % Vérification avant appel
    fprintf(' Appel de trace_erreur2 avec : [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', ...
            k1F, k1B, k2F, k2B, k3F, k3B, k4);

    % Appel de trace_erreur2
    [J1opt, J_rel] = trace_erreur2(E0, texp, xexp, tsw, k1F, k1B, k2F, k2B, k3F, k3B, k4);

    % Gestion des erreurs
    if exist('J1opt', 'var') && ~isempty(J1opt)
        error_stat = J1opt;
        error_stat_rel = J_rel;
    else
        error_stat = inf;
        fprintf(' J1opt non défini ou vide. Erreur assignée à inf.\n');
    end  

    % Calcul des erreurs min-max
    [error_min, error_max, error_min_rel, error_max_rel] = objectif_minmax(E0, texp, xexp, k1F, k1B, k2F, k2B, k3F, k3B, k4);
end


% à faire ce soir : mettre au propre mes avancées
% à faire demain : changer de solution initiale pour voir si on est vrmt
% bloqué dans un min local
% à faire demain : voir ce les conseils par rapport au script Python

