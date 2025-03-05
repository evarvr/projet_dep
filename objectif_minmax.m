function [J_min, J_max,J_min_rel,J_max_rel] = objectif_minmax(E0, texp, xexp, k1F, k1B, k2F, k2B, k3F, k3B, k4)
    % GLOBALS nécessaires pour accéder aux données
    %global texp xexp
    global tref xref
    %global tsw Itsw

    %------------------------------------------
    % Simulation
    %------------------------------------------
    p=[k1F k1B k2F k2B k3F k3B k4];
    [t, x] = simulate1(p, tref);

    %------------------------------------------
    % Calcul de l'erreur min-max
    %------------------------------------------
    J_min = zeros(4, 1);
    J_max = zeros(4, 1);
    rel_error_min = zeros(4, 1);  % Initialisation pour l'erreur relative minimale
    rel_error_max = zeros(4, 1);
    IColonne = [1 2 4 5];
    for i = 1:4
        % Sélection de la colonne à considérer
        Ix = IColonne(i);
        % Sélection des lignes valides
        Iuse = find(~isnan(xref(:, Ix)));
        % Calcul des écarts
        err = abs(xref(Iuse, Ix) - x(Iuse, Ix));
        % Calcul de l'erreur minimale et maximale
        J_min(i) = min(err);
        J_max(i) = max(err);


        rel_err = (xref(Iuse, Ix) - x(Iuse, Ix)) ./ (xref(Iuse, Ix) .* (xref(Iuse, Ix) >= 1e-5) + (xref(Iuse, Ix) < 1.0e-5));
       
        % Calcul de l'erreur relative minimale et maximale
        rel_error_min(i) = min(abs(rel_err));
        rel_error_max(i) = max(abs(rel_err));

    end

    % Pondération des termes
    J_min = 0.25 * J_min(1) + 0.25 * J_min(2) + 0.25 * J_min(3) + 0.25 * J_min(4);
    J_max = 0.25 * J_max(1) + 0.25 * J_max(2) + 0.25 * J_max(3) + 0.25 * J_max(4);
    J_min_rel=0.25 * rel_error_min(1) + 0.25 * rel_error_min(2) + 0.25 * rel_error_min(3) + 0.25 * rel_error_min(4);
    J_max_rel=0.25 * rel_error_max(1) + 0.25 * rel_error_max(2) + 0.25 * rel_error_max(3) + 0.25 * rel_error_max(4);
    J_min_rel = J_min_rel ;
    J_max_rel = J_max_rel ;


end
