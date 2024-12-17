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
k1F_values = linspace(1000, 100000, 3);  
k1B_values = linspace(1000, 10000, 3); 
k2F_values = linspace(0.1,10, 3);
k2B_values = linspace(100,1000, 3);
k3F_values = linspace(1,1000, 3);
k3B_values = linspace(10,1000, 3);
k4_values  = linspace(1,1000, 3);

% Initialisation des variables pour stocker les résultats
best_params = [NaN, NaN]; % Initialisation avec des valeurs NaN
best_error = inf;  % Erreur initiale très élevée

% Boucle sur les combinaisons des valeurs de k1F et k1B
iteration = 0;
for k1F = k1F_values(:)'
    for k1B = k1B_values(:)'
        for k2F=k2F_values(:)'
            for k2B=k2B_values(:)'
                for k3F=k3F_values(:)'
                    for k3B=k3B_values(:)'
                        for k4=k4_values(:)'
                            error = simulate_and_evaluate(E0, texp, xexp, k1F, k1B, k2F, k2B, k3F, k3B, k4);
                            %fprintf('k1F: %.4f, k1B: %.4f, Erreur: %.4f\n', k1F, k1B, error);
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
    %global texp xexp tref xref E0 Es0 Ed0 tsw Itsw xsw kd modele k1F k1B k2F k2B k3F k3B k4 numfich
    

    % Appel du simulateur (en tant que script)
    run('Simulateur.m');  
    
    % Vérification de l'existence des données simulées
    if exist('x1', 'var') && ~isempty(x1)
        % Transposer si nécessaire pour avoir les lignes comme dimension temporelle
        if size(x1, 2) > size(x1, 1)
            x1 = x1';
        end

        % Création d'un vecteur de temps pour les données simulées
        t_sim = linspace(min(texp), max(texp), size(x1, 1));

        % Interpoler les données simulées aux points temporels des données expérimentales
        x1_interp = interp1(t_sim, x1, texp, 'linear', 'extrap');

        % Réduction pour correspondre à xexp en utilisant les 5 premières colonnes
        x1_interp = x1_interp(:, 1:size(xexp, 2));

        % Calcul de l'erreur si les dimensions sont compatibles
        if size(xexp, 2) == size(x1_interp, 2)
            error = sum((xexp - x1_interp).^2, 'all');  % Erreur totale sur tous les points expérimentaux
        else
            error = inf;  % Dimensions incompatibles, attribuer une grande erreur
            fprintf('Incompatibilité de colonnes entre xexp et x1_interp, erreur assignée à inf.\n');
        end
    else
        error = inf;  % Si x1 est vide ou non défini, attribuer une grande erreur
        fprintf('x1 n''existe pas ou est vide, erreur assignée à inf.\n');
    end
end




