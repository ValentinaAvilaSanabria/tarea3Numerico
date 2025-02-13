
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

        % Aplicar valores de bordes solo cuando es necesario
        if i == 1, f_interior(index) = f_interior(index) - bordes2(1, j_img); end
        if i == n, f_interior(index) = f_interior(index) - bordes2(N, j_img); end
        if j == 1, f_interior(index) = f_interior(index) - bordes2(i_img, 1); end
        if j == n, f_interior(index) = f_interior(index) - bordes2(i_img, N); end
    end
end


max_iter = 50; % Menos iteraciones porque Gauss-Seidel converge más rápido
tol = 1e-3; % Margen de error más amplio para reducir tiempo de cómputo

% Inicialización en ceros
u_interior = zeros(n2, 1);


for k = 1:max_iter
    u_prev = u_interior; % Guardamos la iteración anterior para comparar convergencia
    
    for i = 1:n2
        % Actualizar el valor inmediatamente usando los valores ya calculados
        suma = 0;
        if i > 1
            suma = suma + A(i, i-1) * u_interior(i-1); % Valor izquierdo (ya actualizado)
        end
        if i < n2
            suma = suma + A(i, i+1) * u_interior(i+1); % Valor derecho (de iteración anterior)
        end
        if i > n
            suma = suma + A(i, i-n) * u_interior(i-n); % Valor arriba (ya actualizado)
        end
        if i <= n2 - n
            suma = suma + A(i, i+n) * u_interior(i+n); % Valor abajo (de iteración anterior)
        end
        
        u_interior(i) = (f_interior(i) - suma) / A(i, i); % Nueva actualización
    end

    % Criterio de convergencia con margen de error más amplio
    if norm(u_interior - u_prev, inf) < tol
        fprintf(' Convergencia alcanzada en %d iteraciones.\n', k);
        break;
    end
end

if k == max_iter
    fprintf(' El método de Gauss-Seidel no convergió en %d iteraciones.\n', max_iter);
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
title('Imagen reconstruida ');

