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
% INPUTS:
%   A = matriz de coeficientes (matriz de restricciones)
%   b = vector columna de recursos (término independiente restricciones)
%   c = vector fila de costes (coeficientes función objetivo)
% OUTPUTS:
%   x = vector fila solución del problema de programación linea
%   z = escalar valor de la función objetivo para la solución obtenida

    % Determinamos el tamaño de la matriz de restricciones
    [n, m] = size(A);

    % Se añaden variables de holgura con coeficientes igual a cero
    c = horzcat(c, zeros(1, n));

    % Guardamos los índices de las variables que pertenecen a la base
    % Inicialmente las variables de holgura estarán en la base
    base = 1+m:1+m+n;

    % Con las variables de holgura el sistema Ax<=b pasa a ser Ax=b
    % Se procede a definir la tabla de Charles, Cooper y Henderson
    % Esta tabla contendrá en su primera columna los valores de las variables de la base
    aij = horzcat(b, A, eye(n));

    % Definimos las últimas filas de la tabla que determinarán la columna pivote
    % La primera posición de zj contendrá el valor de la función objetivo
    zj = zeros(1, 1+m+n);
    cj_zj = c;
    
    disp(aij)
    cont = 0;
    if all(cj_zj<=0)
        fprintf("ERROR: El problema planteado no tiene solución acotada\n")
        return
    end

    % Mientras que cj-zj sea mayor que cero se repite el algoritmo
    while any(cj_zj > 0)
        % Para obtener la columna pivote k, se busca el valor más grande positivo
        [~, k] = max(cj_zj);
        k = k + 1; % MATLAB indexa a partir de 1

        % Se contempla el caso de dividir por cero
        % MATLAB devolverá infinito pero buscamos el mínimo valor positivo
        a = aij(:, 1)./aij(:, k);
        [~, h] = min(a(a>0));
        fprintf("Pivote(%d, %d) = %.2f\n", h, k, aij(h, k))
        % GAUSS
        for i=1:n
            for j = 1:1+m+n
                if i==h || j==k
                    continue
                end
                aij(i, j) = aij(i, j) - aij(h, j) * aij(i, k) / aij(h, k);
            end
        end
        % Finalmente se modifican la fila y columna pivote
        aij(h, :) = aij(h, :) ./ aij(h, k);
        % Se convierten a cero los elementos en la misma columna que el pivote
        % Y el pivote pasa a ser 1
        aij(:, k) = 0;
        aij(h, k) = 1;



        disp(aij)
        cont = cont+1;
        if cont > 0
            break
        end
    end
end