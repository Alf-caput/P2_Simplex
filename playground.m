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
%   x = vector fila solución del problema de programación lineal
%   z = valor de la función objetivo para la solución obtenida

    % Determinamos el tamaño de la matriz de restricciones
    [n, m] = size(A);

    % A la función objetivo se añaden variables de holgura con coeficientes igual a cero
    cj = horzcat(c, zeros(1, n));

    % Guardamos los índices de las variables que pertenecen a la base
    % Inicialmente las variables de holgura estarán en la base
    base = m+1:m+n; % Inicialmente m variables y a estas se añaden n de holgura

    % Con las variables de holgura el sistema Ax<=b pasa a ser Ax=b
    % Se procede a definir la tabla de Charles, Cooper y Henderson
    % Esta tabla contendrá en su primera columna los valores de las variables de la base
    aij = horzcat(b, A, eye(n));

    % Definimos las últimas filas de la tabla que determinarán la columna pivote
    % La primera posición de zj contendrá el valor de la función objetivo
    % zj se obtiene como el producto escalar aj*cj' 
    % Con cj' coeficientes de las variables en la base
    zj = zeros(1, 1+m+n);
    % NOTA: Los valores de las variables que no están en la base son cero
    % Por eso hacemos aj*cj' en vez de todos los valores por sus respectivos costes

    % cj-zj inicialmente siempre igual a c al ser zj vector de ceros inicialmente
    cj_zj = cj;
    
    if all(cj_zj<=0)
        fprintf("ERROR: El problema planteado no tiene solución acotada\n")
        return
    end

    fprintf("Base inicial:")
    disp(base)
    disp(aij)
    disp(zj)
    disp(cj_zj)
    % Mientras que todos los cj-zj sean mayor que cero se repite el algoritmo
    while any(cj_zj > 0)
        % Obtención columna pivote k
        % Se busca el valor más grande positivo de cj-zj
        [~, k] = max(cj_zj);
        k = k + 1; % MATLAB indexa a partir de 1

        % Obtención fila pivote h
        % Se busca el menor de dividir el valor de las variables que están en la base
        % entre las coordenadas del vector entrará en la base
        a = aij(:, 1)./aij(:, k);
        [~, h] = min(a(a>0));
        % NOTA: Se contempla el caso de dividir por cero
        % MATLAB devolverá infinito pero buscamos el mínimo valor positivo evitando el error

        % Entrará en la base la variable xk
        % Saldrá la variable con posición h en la base
        base(h) = k-1;

        fprintf("Pivote(%d, %d) = %.2f\n", h, k, aij(h, k))
        fprintf("La nueva base es:")
        disp(base)
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
        % Se divide la fila pivote entre el pivote
        aij(h, :) = aij(h, :) ./ aij(h, k);
        % Se hacen cero los elementos en la columna pivote
        aij(:, k) = 0;
        % Excepto el pivote que vale 1
        aij(h, k) = 1;

        % Obtenemos cada zj haciendo el producto escalar aj*cj'
        % Con cj' los coeficientes de las variables en la base
        for j=1:1+m+n
            zj(j) = cj(base)*aij(:, j);
        end
        % Finalmente se obtiene cj-zj
        % Que servirá para calcular el siguiente pivote o detener el algoritmo
        cj_zj = cj-zj(2:end);
        % NOTA: El primer valor de zj no interviene para elegir la columna pivote
        % El primer valor de zj es el valor de la función objetivo para esa base

        % Definimos también una tolerancia
        % En caso de obtener valores muy pequeños los sustituimos por cero
        cj_zj(abs(cj_zj)<=100*eps) = 0;

        disp(aij)
        disp(zj)
        disp(cj_zj)
    end
end