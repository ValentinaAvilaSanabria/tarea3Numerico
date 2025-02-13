
N = size(bordes2,1);
n = N - 2;
n2 = n * n;


e = ones(n2,1);
main_diag  = -4 * e;
side_diag  = 1 * e;
updown_diag = 1 * e;

A = spdiags([updown_diag side_diag main_diag side_diag updown_diag], ...
            [-n -1 0 1 n], n2, n2);
A = sparse(A);


f_interior = bordes2(2:N-1, 2:N-1);
f_interior = f_interior(:);


for i = 1:n
    for j = 1:n
        index = (j-1) * n + i;
        i_img = i + 1;
        j_img = j + 1;

        
        if i == 1, f_interior(index) = f_interior(index) - bordes2(1, j_img); end
        if i == n, f_interior(index) = f_interior(index) - bordes2(N, j_img); end
        if j == 1, f_interior(index) = f_interior(index) - bordes2(i_img, 1); end
        if j == n, f_interior(index) = f_interior(index) - bordes2(i_img, N); end
    end
end


max_iter = 300; % Número máximo de iteraciones
tol = 1e-6; % Tolerancia para la convergencia

% Inicialización en ceros
u_interior = zeros(n2, 1);

D = diag(A); % Extraer la diagonal de A
R = A - diag(D); % Parte restante de A (L + U)


for k = 1:max_iter
    u_new = (f_interior - R * u_interior) ./ D; % Fórmula de Jacobi
    
    % Criterio de convergencia
    if norm(u_new - u_interior, inf) < tol
        fprintf(' Convergencia alcanzada en %d iteraciones.\n', k);
        break;
    end
    
    u_interior = u_new; % Actualizar solución
end

% Verificar si el método no convergió
if k == max_iter
    fprintf(' El método de Jacobi no convergió en %d iteraciones.\n', max_iter);
end


U = bordes2;
U(2:N-1, 2:N-1) = reshape(u_interior, [n, n]);

figure;
surf(U, 'EdgeColor', 'none'); % Gráfica en 3D sin bordes
shading flat; % Suaviza el color sin interpolar
axis ij; % Eje Y en orientación normal
axis equal; % Escalar los ejes uniformemente
view(2); % Vista en 2D desde arriba
colormap bone; % Mapa de colores con más contraste
colorbar;
title('Imagen reconstruida');


figure;
imagesc(U);
axis ij;
axis equal;
colormap gray; % Mapa de colores en escala de grises
colorbar;
title('Imagen reconstruida');