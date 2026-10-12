% System parameters
B = .009;         % Irrigation effectiveness
xref = 0.3;      % Moisture setpoint
xmin = 0.15;      % Minimum moisture
xmax = 0.4;      % Maximum moisture
umax = 10;      % Max irrigation
kE = 0.027;       % Evaporation coefficient

% Time parameters
start_time = datetime('now', 'Format', 'HH'); % Start at current hour
sim_duration_days = 2;                        % 7-day simulation
time_steps_per_hour = 1;                      % 1-hour intervals
total_steps = sim_duration_days * 24;         % Total simulation steps

% MPC parameters
N = 24; % 24-hour prediction horizon
p = 4;  % 4-hour control horizon (6AM-10AM)

% Weather API configuration
apiKey = 'd6bc4ca03f710c49bb0b20df9e686187'; % Replace with actual API key
location = 'Ghaziabad,IN'; % e.g., 'lat=37.7749&lon=-122.4194'

% Initialize simulation
x = zeros(total_steps+1, 1); 
x(1) = 0.4; % Initial moisture
u = zeros(total_steps, 1); 

% Initialize data logging arrays
TA_actual = [];
RH_actual = [];
rain_actual = [];
E_actual = [];

% Simulation loop
for t = 1:total_steps
    current_time = dateshift(start_time, 'start', 'day') + hours(t-1);
    
    % Check irrigation window (6AM-10AM)
    current_hour = hour(current_time);
    irrigation_allowed = (current_hour >= 6) && (current_hour < 10);
    
    % Get weather forecast every 6 hours
    if mod(t-1, 6) == 0 || t == 1
        [TA_forecast, RH_forecast, rain_forecast] = get_weather_forecast(apiKey, location, N);
        
        % Precompute evaporation rates for entire horizon
        E_forecast = kE .* (TA_forecast - (TA_forecast - (100 - RH_forecast)/5)) .* (1 - RH_forecast/100);
    end
    
    % Get current weather conditions
    TA = TA_forecast(1);
    RH = RH_forecast(1);
    rain = rain_forecast(1);
    E = E_forecast(1);
    
    % Log actual weather data
    TA_actual = [TA_actual; TA];
    RH_actual = [RH_actual; RH];
    rain_actual = [rain_actual; rain];
    E_actual = [E_actual; E];
    
    % Simulate sensor reading with noise
    if t == 1
        x_current = x(t); % Initial condition
    else
        x_current = x(t) + 0.01*randn; 
        x_current = max(xmin, min(xmax, x_current)); % Clip to physical limits
    end
    
    % Only calculate control input during allowed hours
    if irrigation_allowed
        % Build prediction matrices
        G = zeros(N, p);
        S = zeros(N, 1);
        D = zeros(N, 1);
        product_terms = 1;
        
        for k = 1:N
            % Update product term for S
            product_terms = product_terms * (1 - E_forecast(k));
            S(k) = product_terms;
            
            % Calculate D(k) with rainfall accumulation
            if k == 1
                D(k) = rain_forecast(k);
            else
                D(k) = sum(rain_forecast(1:k) .* cumprod([1, (1 - E_forecast(1:k-1))]));
            end
            
            % Calculate G matrix
            for j = 1:min(p, k)
                G(k,j) = B * prod(1 - E_forecast(j+1:k));
            end
        end
        
        % Quadratic programming setup
        Q = eye(N);
        H = 2 * G' * Q * G;
        f = 2 * G' * Q * (S*x_current + D - xref);
        
        % Constraints
        A = [eye(p); -eye(p)];
        b = [umax*ones(p,1); zeros(p,1)];
        
        % State constraints (converted to input constraints)
        for i = 1:N
            A = [A; G(i,:)];
            b = [b; xmax - S(i)*x_current - D(i)];
            
            A = [A; -G(i,:)];
            b = [b; -xmin + S(i)*x_current + D(i)];
        end
        
        % Solve QP with fallback
        options = optimoptions('quadprog', 'Display', 'off');
        [u_candidate,~,exitflag] = quadprog(H, f, A, b,[],[],[],[],[],options);
        
        if exitflag <= 0 || isempty(u_candidate)
            u(t) = 0;
            warning('QP failed at t=%d, using u=0', t);
        else
            u(t) = max(0, min(umax, u_candidate(1)));
        end
    else
        u(t) = 0;
    end
    
    % Update state with system dynamics
    if t < total_steps
        x(t+1) = (1 - E)*x_current + B*u(t) + rain;
    end
    
    % Update forecasts (maintain length N)
    E_forecast = [E_forecast(2:end), E_forecast(end)];
    TA_forecast = [TA_forecast(2:end), TA_forecast(end)];
    RH_forecast = [RH_forecast(2:end), RH_forecast(end)];
    rain_forecast = [rain_forecast(2:end), rain_forecast(end)];
end

% Plotting section
figure;
% subplot(4,1,1);
plot(0:total_steps, x, 'b', 0:total_steps, xref*ones(total_steps+1,1), 'r--');
ylabel('Moisture');
legend('Actual', 'Setpoint');

figure;
% subplot(4,1,2);
stairs(0:total_steps-1, u, 'g');
ylabel('Irrigation');
xlabel('Time (hours)');

% % figure;
% subplot(4,1,3);
% plot(0:total_steps-1, TA_actual, 'r', 0:total_steps-1, RH_actual, 'b');
% ylabel('Weather');
% legend('Temperature (C)', 'RH (%)');
% 
% subplot(4,1,4);
% plot(0:total_steps-1, E_actual, 'm');
% ylabel('Evaporation');
% xlabel('Time (hours)');

% Weather API interface function
function [TA, RH, rain] = get_weather_forecast(apiKey, location, N)
    % Simulated API response - replace with actual API call
    t = 0:N-1;
    TA = 20 + 5*sin(2*pi*t/24); % Temperature (C)
    RH = 60 - 20*cos(2*pi*t/24); % Relative humidity (%)
    rain = 0.02 + 0.01*rand(1, N); % Rainfall
end