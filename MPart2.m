%%%%%%%%%%%%%%%%%%%%%%%% Code from Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
%pkg load communications

%% Sampling frequency
global fs
fs = 20000;             % 20KHz
global ts
ts = 1/fs;              % sampling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate stream of random bits (10,000 bits)
N = 10000;                  % number of bits
bits = randi([0 1], 1, N);  % random binary vector

%% Choose a bit rate (bit/s)
Rb = 1000;
Tb = 1/Rb;

%% Sampling frequency
fs = 20000;             % 20KHz
ts = 1/fs;              % sampling period
ds = fs/Rb;             % samples per bit

%% time domain
T = N * Tb;             % total time
t = 0 : ts : (T-ts);    % time vector

y = [];     % amplitude vector

%% Encode the bits in PNRZ
k = 1;
for i = 1:N
    if bits(i) == 1
        for j = 1:ds
            y(k) = 1;
            k = k + 1;
        end
    elseif bits(i) == 0
        for j = 1:ds
            y(k) = -1;
            k = k + 1;
        end
    end
end

figure;
plot(t, y);
grid on;
title('Transmitted Signal', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 14);
ylabel('Amplitude (v)', 'FontSize', 14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Eb_over_No_dB = [-10 -8 -6 -4 -2 0 2 4 6];
BER_s = noise_sweep(t, y, N, ds, bits, Eb_over_No_dB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Theoretical BER %%%%%%%%%%%%%%%%%%%%%%%%%%%

Eb_over_No = [];
BER_t = [];
for i=1:length(Eb_over_No_dB)
    Eb_over_No(i) = 10^(Eb_over_No_dB(i)/10);
    BER_t(i) = qfunc(sqrt(2*Eb_over_No(i)));
end

% Plot BER vs Eb/No for actual & theoretical
figure;
semilogy(Eb_over_No_dB, BER_s, 'b-');
hold on;
semilogy(Eb_over_No_dB, BER_t, 'r-');
grid on;
legend('Actual', 'Theoretical', 'FontSize', 16);
title('BER vs E_b/N_o ', 'FontSize', 16);
xlabel('E_b/N_o (dB)', 'FontSize', 14);
ylabel('BER', 'FontSize', 14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Polar NRZ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume input that the Amplitude of the signal is:
% 1 --> 1 & -1 --> 0, and period T = 1s

PNRZ_sig = {sample([1;1]), sample([1;-1])};
T = PNRZ_sig{1}(1);
c = basis(PNRZ_sig, T);
symbol_energy = constellation_diagram(c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Sampling Function %%%%%%%%%%%%%%%%%%%%%%%%%

function sampled_vector = sample(input_vector)
    global ts;
    y = [];                         % vector of amplitudes at any t
    n = size(input_vector, 2);      % n represents the number of
    t_start = 0;                    % different amplitudes in the pulse
    k = 1;                          % loop iterator
    for i=1:n
        t_end = input_vector(1,i);
        for j=t_start: ts : t_end - ts
            y(k) = input_vector(2,i);
            k = k + 1;
        end
        t_start = t_end;
    end
    sampled_vector = y;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Energy Calc. %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function energy = CalcEnergy(x, T)
    n = length(x);
    delta = 1/n;
    eg = 0;
    for i = 1:n
          eg = eg + x(i)^2 * delta;
    end
    energy = eg * T;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Signal Component %%%%%%%%%%%%%%%%%%%%%%%%%%

function coefficient = CalcCoefficient(s, phi, T)
    n = length(s);
    delta = 1/n;
    coef = 0;
    for i=1:n
        coef = coef + s(i) * phi(i) * delta;
    end
    coefficient = coef * T;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Basis Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%

function sig_comp = basis(input_signals, T)
    global fs;
    global ts;
    len = T * fs;
    t = 0 : ts : (T-ts);    % time vector
    M = length(input_signals);
    basis_vectors = {};      % cell array of basis functions
    phi = [];      % Initialization of phi amplitude vector
    c = [];          % Initialization of coefficients vector

    %% loop through signals
    for num=1:M
        flag = 0;                % flag to check wether g = 0 or not
        N = length(basis_vectors);
        y = input_signals{num};  % amplitude vector of signal

        %% calculating basis function
        if num == 1                 % num represents the signal number
            energy = CalcEnergy(input_signals{1}, T);
            for i = 1:len
                phi(i) = y(i) / sqrt(energy);
            end
            basis_vectors{1} = phi;
            c(1,1) = CalcCoefficient(input_signals{1}, basis_vectors{1}, T);
        else
            sig = input_signals{num};
            g = sig;
            for i=1:N
                phi = basis_vectors{i};
                comp = CalcCoefficient(sig, phi, T);
                if abs(comp) < 1e-10
                    comp = 0;
                end
                c(num, i) = comp;
                g = g - comp * phi;
            end

            for i=1:len
                if abs(g(i)) < 1e-10    % elimnate any errors caused
                    g(i) = 0;         % from the subtraction
                end                   % of floating numbers
                if g(i) ~= 0
                    flag = 1;
                end
            end

            if flag == 1                % if g != 0
                energy = CalcEnergy(g, T);
                for i = 1:len
                    phi(i) = g(i) / sqrt(energy);
                end
                basis_vectors{N+1} = phi;
                c(num,num) = sqrt(energy);
            end
        end
    end
    sig_comp = c;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Constellation Diagram %%%%%%%%%%%%%%%%%%%%%%%%

function symbol_energy = constellation_diagram(matrix)

    % Extract coordinates
    [r,c] = size(matrix);
    x = matrix(:, 1);

    % Check if 3d or 2d or 1d
    if c == 3
        y = matrix(:, 2);
        z = matrix(:, 3);
        energy = x.^2 + y.^2 + z.^2;

        figure;
        plot3(x, y, z, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        grid on;
        title('Constellation Diagram', 'FontSize', 16);
        xlabel('Ф_1', 'FontSize', 14);
        ylabel('Ф_2', 'FontSize', 14);
        zlabel('Ф_3', 'FontSize', 14);

        % add labels for each point
        for i = 1:r
            text(x(i), y(i), z(i), sprintf('   S_%d', i),...
                 'FontSize', 14);
        end
    elseif c == 2
        y = matrix(:, 2);
        energy = x.^2 + y.^2;

        % plot
        figure;
        plot(x, y, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        grid on;
        title('Constellation Diagram', 'FontSize', 16);
        xlabel('Ф_1', 'FontSize', 14);
        ylabel('Ф_2', 'FontSize', 14);

        % add labels for each point
        for i = 1:r
            text(x(i), y(i), sprintf('   S_%d', i), 'FontSize', 14);
        end
    else
        energy = x.^2;
        figure;
        plot(x, zeros(size(x)), 'o', 'MarkerSize', 8,...
            'MarkerFaceColor', 'b');
        grid on;
        title('Constellation Diagram', 'FontSize', 16);
        xlabel('Ф_1', 'FontSize', 14);

        % add labels for each point
        for i = 1:r
            text(x(i), 0, sprintf('   S_%d', i), 'FontSize', 14);
        end
    end

    % Return the vector containing the energy of each symbol
    symbol_energy = energy;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%% BER Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = calculate_BER(n,txb,rxb)
  error_count = 0;
  for i = 1:n
    if rxb(i) ~= txb(i)
      error_count = error_count + 1;
    end
  end
  val = error_count / n;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% Decision Device %%%%%%%%%%%%%%%%%%%%%%%%%%

function rxb = decision_device(n, samples_per_bit, amplitude_vector)
    k = samples_per_bit;
    for i = 1:n
        v = amplitude_vector(k);
        k = k + samples_per_bit; % sample every Tb
        if real(v) < 0
          rxb(i) = 0;
        else
          rxb(i) = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Noise Sweep %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BER_s = noise_sweep(time_vector, a_vector, n,...
                            samples_per_bit, txb, Eb_over_No)
    len = length(Eb_over_No);
    snr = Eb_over_No + 10 * log10(n);
    c_v = a_vector + 0.0000000000000000000000000000000000000000000001i;
    pn = zeros(2,1);
    pn(1,1) = 1;
    pn(2,1) = -1;
    
    for i = 1:len
        %% Signal Power = (N * amplitude ^ 2 * Tb) / Tb = N
        s = awgn(c_v, snr(i), 10 * log10(n));
        rxb = decision_device(n, samples_per_bit, s);
        BER_s(i) = calculate_BER(n, txb, rxb);

%%%%%%%%%%%%%%%%%% Complex Constellation Diagram %%%%%%%%%%%%%%%%%%%%%

        swcn = s(samples_per_bit : samples_per_bit : n*samples_per_bit);
        [r,c] = size(pn);
        x = pn(:, 1);
        energy = x.^2;
        figure;
        scatter(real(swcn), imag(swcn),7,...
            'd','MarkerFaceColor','r','MarkerEdgeColor','r');
        xlim([-5 5]);
        ylim([-5 5]);
        hold on;
        plot(x, zeros(size(x)), 'o', 'MarkerSize', 8,...
            'MarkerFaceColor', 'b');
        grid on;
        title(['Constellation Diagram with Complex Noise at E_b/N_o = ' ... 
        num2str(Eb_over_No(i)) ' dB'], 'FontSize', 16);
        xlabel('Real', 'FontSize', 14);
        ylabel('Imaginary', 'FontSize', 14);
        % add labels for each point
        for ord = 1:r
            text(x(ord), 0, sprintf('   S_%d', ord), 'FontSize', 14);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end


    % Plot noisy signal
    figure;
    plot(time_vector, s);
    grid on;
    title('Received Signal', 'FontSize', 16);
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Amplitude (v)', 'FontSize', 14);

    % Plot BER vs Eb/No for simple decision device
    figure;
    semilogy(Eb_over_No, BER_s);
    grid on;
    title('BER vs E_b/N_o', 'FontSize', 16);
    xlabel('E_b/N_o (dB)', 'FontSize', 14);
    ylabel('BER', 'FontSize', 14);

end
