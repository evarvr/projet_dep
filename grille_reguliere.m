%clear all;
close all;

tic; % debut du temps

global texp xexp tref xref E0 Es0 Ed0 tsw Itsw xsw kd modele k1F k1B k2F k2B k3F k3B k4 numfich


% Initialisation du programme
modele = 'bbo'; % Choix du modèle
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
best_params = [NaN, NaN]; % Initialisation avec des valeurs NaN
best_error = inf;  % Erreur initiale très élevée

% Boucle sur les combinaisons des valeurs de k1F et k1B
iteration = 0;
error_list = []; % Initialiser la liste des erreurs
for k1F = k1F_values(:)'
    for k1B = k1B_values(:)'
        for k2F=k2F_values(:)'
            for k2B=k2B_values(:)'
                for k3F=k3F_values(:)'
                    for k3B=k3B_values(:)'
                        for k4=k4_values(:)'
                            error = simulate_and_evaluate(E0, texp, xexp, k1F, k1B, k2F, k2B, k3F, k3B, k4);
                            error_list = [error_list; error];
                            iteration = iteration + 1;
                            if error < best_error
                                best_error = error;
                                best_params = [k1F, k1B, k2F,k2B,k3F,k3B,k4];
                            end
                        end
                    end
                end
            end 
        end
    end
end


% Vérifier si une combinaison valide a été trouvée
if ~isnan(best_params(1))
    fprintf('Nombre d''itérations : %d\n', iteration);
    fprintf('Meilleure combinaison de paramètres : k1F = %.4f, k1B = %.4f, k2F=%.4f, k2B=%.4f,k3F=%.4f, k3B=%.4f, k4=%.4f, Erreur = %.4f\n', best_params(1), best_params(2), best_params(3),best_params(4),best_params(5), best_params(6),best_params(7),best_error);
else
    fprintf('Aucune combinaison de paramètres n''a donné une erreur finie.\n');
end

toc; % fin de la mesure du temps


% Fonction de simulation et évaluation
function error = simulate_and_evaluate(E0, texp, xexp, k1F, k1B, k2F, k2B, k3F, k3B, k4)

    % Appel du simulateur (Simulateur.m)
    run('Simulateur.m');

    % Si le simulateur renvoie J1opt directement, utiliser cette valeur comme erreur
    if exist('J1opt', 'var') && ~isempty(J1opt)
        error = J1opt; % Utiliser l'erreur globale calculée par le simulateur
    else
        % Si J1opt n'est pas défini, attribuer une erreur infinie
        error = inf;
        fprintf('J1opt n''est pas défini ou est vide. Erreur assignée à inf.\n');
    end
end


% Calculer la somme des données expérimentales pour normalisation
sum_exp = sum(xexp(:));

% Calculer l'erreur relative
relative_error = best_error / sum_exp;

% Afficher l'erreur relative pour mieux comprendre
fprintf('Erreur relative : %.4f (ou %.2f%%)\n', relative_error, relative_error * 100);

%figure;
%plot(1:length(error_list), error_list, '-o', 'LineWidth', 1.5);
%xlabel('Itérations');
%ylabel('Erreur');
%title('Évolution de l''erreur au fil des itérations');
%grid on;


