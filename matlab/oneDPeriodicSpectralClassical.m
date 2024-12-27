% Parameters
time = 2;

C = 4;
D = 1;
S = -0.2;
wx = 2;
L = wx * pi;

% Number of qubits for physical space
nq = 8;
nx = 2^nq;

% Spatial grid
x = linspace(-L/2, L/2, nx);

% Wavenumbers used in initial condition
k_s = [1, 3];
k_c = 2;

u_init = sum(sin(k_s(:) * x), 1) + sum(cos(k_c(:) * x), 1);

% Fourier domain setup
kx = fftshift((0 : nx-1) - nx/2) * (2 * pi / L);
fx = fft(u_init);

rhs = @(t, y) reshape((-1i*C*kx(:) - D*kx(:).^2 + S).*y, [], 1);

sol = ode23(@(t, y) rhs(t, y), [0, time], fx);
t_vals = sol.x;
y_vals = sol.y;

results = ifft(y_vals, [], 1);

% Plot results
figure;
hold on;

for ii = 1:1000
    idx = 2*(ii - 1) + 1;
    t_curr = t_vals(ii);

    ut_s = sum(sin(k_s(:) .* (x - C * t_curr)) ...
              .* exp( (-D * k_s(:).^2 + S) * t_curr ), 1 );
    
    ut_c = sum( cos(k_c(:) .* (x - C * t_curr)) ...
              .* exp( (-D * k_c(:).^2 + S) * t_curr ), 1 );

    ut = ut_s + ut_c;

    plot(x, ut, '-k', 'DisplayName', 'analytic');
    plot(x, real(results(:, ii)), '-.r', 'DisplayName', 'numeric');
end

legend
xlim([-L/2, L/2]);