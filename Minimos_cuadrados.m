clc; clear; close all;

data = load('Datos.txt'); % Cargar el archivo de texto (2 columnas: x y y)
x = data(:,1); % Extraer valores de x
y = data(:,2); % Extraer valores de y

%% Matriz A
n = length(x);
m = 5; 
A = zeros(n, m+1);

for i = 1:m+1
    A(:, i) = x.^(i-1);
end

%% Ecuaciones Normales
AtA = A.' * A;
AtY = A.' * y;
c_normal = AtA \ AtY; 

%% Factorización QR
[Q, R] = qr(A); 
c_qr = R \ (Q' * y); 

%% Verdaderos
c_certificados = ones(m+1, 1); 

% Residuo
residual_normal = norm(A * c_normal - y, 2);
residual_qr = norm(A * c_qr - y, 2);

% Errores relativos
error_relativo_normal = abs((c_normal - c_certificados) ./ c_certificados);
error_relativo_qr = abs((c_qr - c_certificados) ./ c_certificados);

%% Resultados
disp('Coeficientes obtenidos (Ecuaciones Normales):');
disp(c_normal);
disp('Coeficientes obtenidos (QR):');
disp(c_qr);
disp('Coeficientes certificados:');
disp(c_certificados);

disp('Residual ||Ac - y||_2 (Ecuaciones Normales):');
disp(residual_normal);
disp('Residual ||Ac - y||_2 (QR):');
disp(residual_qr);

disp('Diferencia relativa (Ecuaciones Normales):');
disp(error_relativo_normal);
disp('Diferencia relativa (QR):');
disp(error_relativo_qr);

%% Graficas
x_fit = linspace(min(x), max(x), 100); % Puntos para la curva ajustada
y_fit_normal = polyval(flip(c_normal'), x_fit); % Evaluar polinomio con ecuaciones normales
y_fit_qr = polyval(flip(c_qr'), x_fit); % Evaluar polinomio con QR

figure;
scatter(x, y, 'bo', 'filled'); hold on; % Datos originales
plot(x_fit, y_fit_normal, 'r-', 'LineWidth', 2); % Ajuste con ecuaciones normales
plot(x_fit, y_fit_qr, 'g--', 'LineWidth', 2); % Ajuste con QR
legend('Datos', 'Ajuste (Ecuaciones Normales)', 'Ajuste (QR)');
xlabel('x');
ylabel('y');
title('Ajuste Polinomial de Grado 5');
grid on;