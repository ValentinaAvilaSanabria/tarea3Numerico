clc; clear; close all;

%% Matriz A y el vector b
A = [  3  -1  -1   0   0;
      -1   4   0  -2   0;
      -1   0   3  -1   0;
       0  -2  -1   5  -1;
       0   0   0  -1   2];

b = [  2;
     -26;
      3;
     47;
    -10];

n = length(b); 

%% Matriz de los métodos de Jacobi, Gauss-Seidel y Sobrerrelajación
D = diag(diag(A)); % Matriz diagonal de A
L = tril(A, -1);   % Parte triangular inferior de A
U = triu(A,  1);   % Parte triangular superior de A

J = -D \ (L + U); % Matriz de Jacobi
S = -(D + L) \ U; % Matriz de Gauss-Seidel

rho_J = max(abs(eig(J))); % Radio espectral de Jacobi
rho_S = max(abs(eig(S))); % Radio espectral de Gauss-Seidel

fprintf('Radio espectral de Jacobi: %.4f\n', rho_J);
fprintf('Radio espectral de Gauss-Seidel: %.4f\n', rho_S);

%% parámetro óptimo de sobrerrelajación ω*
omega_opt = 2 / (1 + sqrt(1 - rho_S^2)); 
fprintf('Parámetro óptimo de sobrerrelajación ω*: %.2f\n', omega_opt);

%% Cantidad de iteraciones
tol = 1e-6; % Tolerancia
max_iter = 1000;

% Método de Gauss-Seidel
[x_gs, iter_gs] = gauss_seidel(A, b, tol, max_iter);

% Método de Sobrerrelajación con ω*
[x_sor, iter_sor] = sor(A, b, omega_opt, tol, max_iter);

fprintf('Iteraciones necesarias en Gauss-Seidel: %d\n', iter_gs);
fprintf('Iteraciones necesarias en Sobrerrelajación (ω*): %d\n', iter_sor);

iter_diff = iter_gs - iter_sor;
fprintf('Reducción de iteraciones con Sobrerrelajación: %d menos que Gauss-Seidel\n', iter_diff);

%% Mejora en tolerancia
tol = 1e-7;
% Método de Gauss-Seidel
[x_gs, iter_gs] = gauss_seidel(A, b, tol, max_iter);

% Método de Sobrerrelajación con ω*
[x_sor, iter_sor] = sor(A, b, omega_opt, tol, max_iter);

fprintf('Iteraciones necesarias para mejorar en un decimal en Gauss-Seidel: %d\n', iter_gs);
fprintf('Iteraciones necesarias para mejorar en un decimal en Sobrerrelajación (ω*): %d\n', iter_sor);


%% Función para el método de Gauss-Seidel
function [x, iter] = gauss_seidel(A, b, tol, max_iter)
    n = length(b);
    x = zeros(n,1); % Solución inicial
    iter = 0;
    
    for k = 1:max_iter
        x_old = x;
        for i = 1:n
            sum1 = A(i,1:i-1) * x(1:i-1);
            sum2 = A(i,i+1:n) * x_old(i+1:n);
            x(i) = (b(i) - sum1 - sum2) / A(i,i);
        end
        iter = iter + 1;
        if norm(x - x_old, inf) < tol
            break;
        end
    end
end

%% Función para el método de Sobrerrelajación (SOR)
function [x, iter] = sor(A, b, omega, tol, max_iter)
    n = length(b);
    x = zeros(n,1); % Solución inicial
    iter = 0;
    
    for k = 1:max_iter
        x_old = x;
        for i = 1:n
            sum1 = A(i,1:i-1) * x(1:i-1);
            sum2 = A(i,i+1:n) * x_old(i+1:n);
            x(i) = (1 - omega) * x_old(i) + (omega * (b(i) - sum1 - sum2) / A(i,i));
        end
        iter = iter + 1;
        if norm(x - x_old, inf) < tol
            break;
        end
    end
end