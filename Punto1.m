% Limpiar variables y pantalla
clear; clc;

% Matriz dada
A = [0 -20 14; -3 27 4; -4 11 2];

% Dimensiones de A
[m, n] = size(A);

%% **Gram-Schmidt**
Q_gs = zeros(m, n);
R_gs = zeros(n, n);

for j = 1:n
    v = A(:, j);
    for i = 1:j-1
        R_gs(i, j) = Q_gs(:, i)' * A(:, j);
        v = v - R_gs(i, j) * Q_gs(:, i);
    end
    R_gs(j, j) = norm(v);
    if R_gs(j, j) ~= 0
        Q_gs(:, j) = v / R_gs(j, j);
    end
end

%% **Householder**
R_h = A;
Q_h = eye(m);

for k = 1:n
    % Seleccionar el subvector desde la fila k en adelante
    x = R_h(k:m, k);
    
    % Evitar errores numéricos con signo correcto
    sigma = norm(x);
    if x(1) >= 0
        sigma = -sigma;
    end
    
    % Crear vector de reflexión v
    v = x;
    v(1) = v(1) - sigma;
    
    % Evitar divisiones por cero
    v_norm = norm(v);
    if v_norm ~= 0
        v = v / v_norm;
    end
    
    % Construcción de la matriz de reflexión H_k
    Hk = eye(m);
    Hk(k:m, k:m) = eye(length(v)) - 2 * (v * v');
    
    % Aplicar transformación de Householder
    R_h = Hk * R_h;
    Q_h = Q_h * Hk';
end

% Asegurar que R_h sea triangular superior
R_h = triu(R_h);

%% **Mostrar resultados**
disp('Q obtenida por Gram-Schmidt:');
disp(Q_gs);
disp('R obtenida por Gram-Schmidt:');
disp(R_gs);

disp('Q obtenida por Householder:');
disp(Q_h);
disp('R obtenida por Householder:');
disp(R_h);


