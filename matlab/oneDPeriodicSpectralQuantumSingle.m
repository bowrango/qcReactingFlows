% Parameters
time = 0;

C = 4;
D = 1;
S = -0.2;
wx = 2;
L = wx * pi;

% Number of qubits for physical space
nq = 4;
nx = 2^nq;

% dimension of the auxillary variable
nqy = 8;
ny = 2^nqy;

wy = 8;

% Define x
x = linspace(-L/2, L/2, nx); % nx is the number of points

% Define wave numbers
k_s = [1, 3]; % Sine components
k_c = 2;    % Cosine components

% Compute sine and cosine components and sum them
u_init = sum(sin(k_s' * x), 1) + sum(cos(k_c' * x), 1);

% Normalize u_init
u_init_norm = sqrt(sum(u_init .^ 2)); % Compute the norm
u0 = u_init / u_init_norm; % Normalize to unit vector

% auxillary variable

shift = 0;
alpha = 10;
bzone = 0;

% Define y
y = linspace(-pi * (wy/2 + shift), pi * (wy/2 - shift), ny); % ny is the number of points

% Define bottom value
bottom = exp(-pi * (wy/2 - shift));

% Compute v_init
v_init = zeros(1, ny);

% Apply conditions
cond1 = (y < -pi * (wy/2 + shift - bzone)) | (y > pi * (wy/2 - shift - bzone));
cond2 = (y >= -pi * (wy/2 + shift - bzone)) & (y < -pi * shift);
cond3 = (y >= -pi * shift);
v_init(cond1) = 0;
v_init(cond2) = max(exp(alpha * (y(cond2) + shift * (1 + 1/alpha) * pi)), bottom);
v_init(cond3) = exp(-y(cond3));

% Perform FFT on v_init
fv = fft(v_init, ny) / sqrt(ny);

% Compute frequency components
kv = fftshift((0:ny-1)/ny - 0.5) * ny * 2 / wy;

% Fourier domain setup
kx = fftshift((0 : nx-1) - nx/2) * (2 * pi / L);

out = zeros(2^nq, 1, 'like', 1j);
for jj = 1:length(kv)

    % wavenumbers
    H = diag(C*kx+D*kv(jj)*((kx).^2)-S*kv(jj));
    U = expm(1j*H*time);

    gates = [
        initGate(1:nq, u0*fv(jj)/abs(fv(jj)))
        qftGate(1:nq)
        unitaryGate(1:nq, U)
        inv(qftGate(1:nq))
        ];
    qc = quantumCircuit(gates);

    s = simulate(qc).Amplitudes;

    out = out +  s.*abs(fv(jj)) .* exp(1j*pi*kv(jj)*wy/2);
end

out = out * u_init_norm;
out = out * exp(y(round(ny / 2)));

% Compute analytical solution directly
uts = sum(sin(k_s' * (x - C * time)) .* exp((-D * (k_s.^2)' + S) * time), 1);
utc = sum(cos(k_c' * (x - C * time)) .* exp((-D * (k_c.^2)' + S) * time), 1);
ut = uts + utc;

figure;
hold on;
plot(x, ut, '-k', 'DisplayName', 'Analytic');
plot(x, real(out), '-.r', 'DisplayName', 'Quantum');
legend

