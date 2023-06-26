% Размеры цилиндра
R = 0.1; % метры
Z = 0.5; % метры

% Размер и количество пузырьков, расположенных внутри цилиндра
n_bubbles = 20;
r_bubble_max = 0.03; % метры
r_bubble_min = r_bubble_max / 5;
r_bubbles = rand(n_bubbles, 1) * (r_bubble_max - r_bubble_min) + r_bubble_min;
x_bubbles = rand(n_bubbles, 1) * (R - r_bubbles);
y_bubbles = rand(n_bubbles, 1) * Z;
z_bubbles = rand(n_bubbles, 1) * (R - r_bubbles);

% Создание цилиндра
theta = linspace(0, 2*pi, 50);
z = linspace(0, Z, 10);
[theta, z] = meshgrid(theta, z);
x = R * cos(theta);
y = R * sin(theta);
C = ones(size(x));
surf(x,y,z,C, 'FaceAlpha',0.5);
hold on

% Создание пузырьков
for i = 1:n_bubbles
    [bx,by,bz] = sphere(20);
    bx = bx * r_bubbles(i) + x_bubbles(i);
    by = by * r_bubbles(i) + y_bubbles(i);
    bz = bz * r_bubbles(i) + z_bubbles(i);
    surf(bx, by, bz, 'FaceAlpha', 0.8);
end

axis equal
xlabel('X, м');
ylabel('Y, м');
zlabel('Z, м');
title('Цилиндр с пузырьками');
