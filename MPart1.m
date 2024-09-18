%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
%pkg load communications

%% Sampling frequency
global fs 
fs = 20000;             % 20KHz
global ts 
ts = 1/fs;              % sampling 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% User Prompts %%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = input("How many signals do you want to input?\n");
T = input("Input Signal period\n");
t = 0 : ts : (T-ts);    % time vector
len = length(t);
signals_vectors = {};    % cell array of user inputs

% input format is a 2xP vector, where P is the number of different regions
% in the signal. row 1 represents time & row 2 represents amplitude
% To clarify, signal 1 in example 1 will be: [T/3 T; 1 0]
for i=1:M
    signals_vectors{i} = sample(input("Input Signal matrix\n"));
end

%% Plotting of signals
for i=1:M
    figure;
    plot(t, signals_vectors{i});
    grid on;
    title(['S_' num2str(i) ' (t)'], 'FontSize', 16);
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Amplitude (V)', 'FontSize', 14);
end

%%%%%%%%%%%%%%%%%%%%%%% Gram-Schmidt procedure %%%%%%%%%%%%%%%%%%%%%%%
c = basis(signals_vectors, T);
symbol_energy = constellation_diagram(c)
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
            %% Plot basis functions
            figure;
            plot(t, phi);
            grid on;
            title(['Ф_1 (t)'], 'FontSize', 16);
            xlabel('Time (s)', 'FontSize', 14);
            ylabel('Amplitude (V)', 'FontSize', 14);
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
                if abs(g(i)) < 1e-10    % eliminate any errors caused
                    g(i) = 0;         % from the subtraction
                end
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
                %% Plot basis functions
                figure;
                plot(t, phi);
                grid on;
                title(['Ф_' num2str(N+1) ' (t)'], 'FontSize', 16);
                xlabel('Time (s)', 'FontSize', 14);
                ylabel('Amplitude (V)', 'FontSize', 14);
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



