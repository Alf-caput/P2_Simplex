clear, clc
c = [5, 4, 3, 1];
A = [
    1, 1, 1, 1;
    3, 2, 3, 1;
    2, 1, 1, 4
];
b = [
    12;
    5;
    7
];
simplex(c, A, b)

function [x, z] = simplex(c, A, b)
% A = matriz de coeficientes parte izquierda restricciones
% b = vector columna de recursos
% c = vector fila de costes
    [n, m] = size(A);
    % Se añaden variables de holgura con coeficientes igual a cero
    c = horzcat(c, zeros(1, n));
    tabla_CCH = horzcat(b, A, eye(n));
    zj = zeros(1, 1+m+n);
    cj_zj = c;
    
    disp(c)
    disp(tabla_CCH)
    disp(zj)
    disp(cj_zj)
    
    while any(cj_zj > 0)
        [~, j_pivote] = max(cj_zj);
        j_pivote = j_pivote + 1; % MATLAB indexa a partir de 1
        [~, i_pivote] = min(tabla_CCH(:, 1)./tabla_CCH(:, j_pivote));
        pivote = tabla_CCH(i_pivote, j_pivote);
        % GAUSS
        % Se convierten a cero los elementos en la misma columna que el pivote
        % Y el pivote pasa a ser 1
        %tabla_CCH()
        break
    end
    fprintf("Pivote(%d, %d) = %.2f\n", j_pivote, i_pivote, pivote)
end