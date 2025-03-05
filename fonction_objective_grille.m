%clear all;
close all;

tic; % debut du temps

global texp xexp tref xref E0 Es0 Ed0 tsw Itsw xsw kd modele k1F k1B k2F k2B k3F k3B k4 numfich


% Initialisation du programme
modele = 'bbpp'; % Choix du modèle
E0 = 10;         % E0 fixé
Es0 = 0.61 * E0; % Valeur dérivée d'E0
Ed0 = 0.39 * E0; % Valeur dérivée d'E0

% Lecture du numéro de fichier UNE SEULE FOIS
numfich = input('Numero de fichier : ','s');
filepath = ['data/mesures-' numfich '.dat'];
if exist(filepath, 'file') ~= 2
    error('Le fichier %s n''existe pas.', filepath);
end

% Charger les données une seule fois
A = load(filepath);
texp = A(:,1);
xexp = A(:,2:end);
tsw = input(['  Temps de commutation du modele sur [0 ; ' num2str(texp(end)) ' ] : ']);

% Afficher les dimensions pour vérification
disp('Dimensions des données expérimentales :');
disp(size(xexp));

% Définir les plages des paramètres à optimiser pour k1F et k1B uniquement
k1F_values = linspace(0.01, 5000, 4);  
k1B_values = linspace(0.01, 5000, 4); 
k2F_values = linspace(0.01,5000, 4);
k2B_values = linspace(0.01,5000, 4);
k3F_values = linspace(0.01,5000, 4);
k3B_values = linspace(0.01,5000, 4);
k4_values  = linspace(0.01,5000, 4);

% Initialisation des variables pour stocker les résultats
best_params_stat = [NaN, NaN]; % Meilleurs paramètres pour erreur statistique
best_params_min = [NaN, NaN];  % Meilleurs paramètres pour erreur minimale
best_params_max = [NaN, NaN];  % Meilleurs paramètres pour erreur maximale

best_error_stat = inf; % Initialisation pour l'erreur statistique
best_error_min = inf;  % Initialisation pour l'erreur minimale
best_error_max = inf;  % Initialisation pour l'erreur maximale
best_error_stat_rel = inf; 
best_error_min_rel = inf;  
best_error_max_rel = inf;  



iteration = 0; % Compteur d'itérations

% Boucle sur les combinaisons des paramètres
for k1F = k1F_values(:)'
    for k1B = k1B_values(:)'
        for k2F = k2F_values(:)'
            for k2B = k2B_values(:)'
                for k3F = k3F_values(:)'
                    for k3B = k3B_values(:)'
                        for k4 = k4_values(:)'
                            % Calcul de l'erreur pour les paramètres actuels
                            [error_stat, error_min, error_max,error_stat_rel, error_min_rel ,error_max_rel] = simulate_and_evaluate(E0, texp, xexp, k1F, k1B, k2F, k2B, k3F, k3B, k4);

                            % Minimiser l'erreur statistique
                            if error_stat < best_error_stat
                                best_error_stat = error_stat;
                                best_error_stat_rel=error_stat_rel;
                                best_params_stat = [k1F, k1B, k2F, k2B, k3F, k3B, k4];
                            end

                            % Minimiser l'erreur minimale
                            if error_min < best_error_min
                                best_error_min = error_min;
                                best_error_min_rel=error_min_rel;
                                best_params_min = [k1F, k1B, k2F, k2B, k3F, k3B, k4];
                            end

                            % Minimiser l'erreur maximale
                            if error_max < best_error_max
                                best_error_max = error_max;
                                best_error_max_rel=error_max_rel;
                                best_params_max = [k1F, k1B, k2F, k2B, k3F, k3B, k4];
                            end

                            iteration = iteration + 1;
                        end
                    end
                end
            end
        end
    end
end

% Affichage des résultats
fprintf('Nombre total d''itérations : %d\n', iteration);

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



%{

figure;

B1   = x1(:,1);
AR1  = x1(:,2);
R1   = x1(:,3);
AB1  = x1(:,4);
ABA1 = x1(:,5);

% Tracé des données expérimentales
plot(texp, xexp(:,1), 'o-', 'DisplayName', 'B1 exp', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(texp, xexp(:,2), 'o-', 'DisplayName', 'AR1 exp', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(texp, xexp(:,3), 'o-', 'DisplayName', 'R1 exp', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(texp, xexp(:,4), 'o-', 'DisplayName', 'AB1 exp', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(texp, xexp(:,5), 'o-', 'DisplayName', 'ABA1 exp', 'LineWidth', 1.5, 'MarkerSize', 8);

% Tracé des données simulées interpolées
plot(t1, B1, '--', 'DisplayName', 'B1 sim', 'LineWidth', 1.5);
plot(t1, AR1, '--', 'DisplayName', 'AR1 sim', 'LineWidth', 1.5);
plot(t1, R1, '--', 'DisplayName', 'R1 sim', 'LineWidth', 1.5);
plot(t1, AB1, '--', 'DisplayName', 'AB1 sim', 'LineWidth', 1.5);
plot(t1, ABA1, '--', 'DisplayName', 'ABA1 sim', 'LineWidth', 1.5);

% Ajustements de la figure
legend('Location', 'Best', 'FontSize', 12);
title('Comparaison entre données expérimentales et simulées', 'FontSize', 14);
xlabel('Temps (unités)', 'FontSize', 12);
ylabel('Concentration ou mesure (unités)', 'FontSize', 12);
grid on;
hold off;


% Calculer la somme des données expérimentales pour la normalisation
sum_exp = sum(xexp(:));



% Fonction de simulation et évaluation pour trois fonctions objectives

function [error_stat, error_min, error_max] = simulate_and_evaluate(E0, texp, xexp, k1F, k1B, k2F, k2B, k3F, k3B, k4)
    % Appel du simulateur (Simulateur.m)
    run('trace_erreur.m');
    
        if exist('J1opt', 'var') && ~isempty(J1opt)
            error_stat = J1opt; % Utiliser l'erreur quadratique
        else
            error_stat = inf;
            fprintf('J1opt n''est pas défini ou est vide. Erreur assignée à inf.\n');
        end   

        % Appeler une fonction objectif min-max
        [error_min,error_max] = objectif_minmax(E0, texp, xexp, k1F, k1B, k2F, k2B, k3F, k3B, k4);
   
end


run('trace_erreur.m'); 
    
figSP = figure;
subplot(2,1,1)
%figS = figure;
plot(texp,Bexp,'ro')
hold on
plot(t1,B1,'r-')
plot(texp,ARexp,'bv')
plot(t1,AR1,'b--')
%plot(t2,B2,'r-')
%plot(t2,AR2,'b--')
hold off
grid on
axis([0 tsw -Inf Inf])
line([tsw tsw],[0 Inf])
xlabel('t (h)')
ylabel('Substrats (M)')
legend('B^{exp}','B','AR^{exp}','AR')

%Graphique (P)
subplot(2,1,2)
%figP = figure;
plot(texp,ABexp,'ro')
hold on
plot(t1,AB1,'r-')
%plot(t2,AB2,'r-')
plot(texp,ABAexp,'bv')
plot(t1,ABA1,'b--')
plot(t1,R1,'m-.')
hold off
grid on
axis([0 tsw -Inf Inf])
line([tsw tsw],[0 Inf])
xlabel('t (h)')
ylabel('Produits (M)')
legend('AB^{exp}','AB','ABA^{exp}','ABA','R')

figE = figure;
subplot(2,1,1)
plot(t1,E1,'r-')
ylabel('Enzyme total')
grid on

subplot(2,1,2)
if strcmp(modele,'synthese1')
    plot(t1,EA,'r-')
    legend('EA')
elseif strcmp(modele,'synthese2')
    plot(t1,EARB,'r-',t1,EARAB,'b--')
    legend('EARB','EARAB')
elseif strcmp(modele,'synthese_mixte')
    plot(t1,EA,'r-',t1,EARAB,'b--')
    legend('EA','EARAB')
end
xlabel('t (h)')
grid on
ylabel('Complexe Enzyme-Substrat')
axis([0 tsw -Inf Inf])


figR = figure;
subplot(2,1,1)
   plot(t1,r1F,'r-^','MarkerSize',3)
   hold on
   plot(t1,r1B,'r--.')
   plot(t1,r2F,'b-<','MarkerSize',3)
   plot(t1,r2B,'b-.')
   hold off
   grid on
   ylabel('r (M/h)')
   legend('r_{1F}','r_{1B}','r_{2F}','r_{2B}')
   
subplot(2,1,2)
   plot(t1,r3F,'r-v','MarkerSize',3)
   hold on
   plot(t1,r3B,'r--.')
   plot(t1,r4, 'b->','MarkerSize',3)
   grid on
   xlabel('t (h)')
   ylabel('r (M/h)')
   legend('r_{3F}','r_{3B}','r_{4}')
   hold off

 % Identifier la longueur minimale
min_length = min(length(texp), length(t1));

% Tronquer les données si nécessaire
texp = texp(1:min_length);
xexp = xexp(1:min_length, :); % Toutes les colonnes expérimentales
t1 = t1(1:min_length);
B1 = B1(1:min_length);
AR1 = AR1(1:min_length);
R1 = R1(1:min_length);
AB1 = AB1(1:min_length);
ABA1 = ABA1(1:min_length);

% Créer une table alignée
figure2_data = table();
figure2_data.t_exp = texp;
figure2_data.B_exp = xexp(:, 1);
figure2_data.AR_exp = xexp(:, 2);
figure2_data.R_exp = xexp(:, 3);
figure2_data.AB_exp = xexp(:, 4);
figure2_data.ABA_exp = xexp(:, 5);

figure2_data.t_sim = t1;
figure2_data.B_sim = B1;
figure2_data.AR_sim = AR1;
figure2_data.R_sim = R1;
figure2_data.AB_sim = AB1;
figure2_data.ABA_sim = ABA1;

% Afficher la table
disp('Tableau des données utilisées pour la figure 2 :');
disp(figure2_data);


%}