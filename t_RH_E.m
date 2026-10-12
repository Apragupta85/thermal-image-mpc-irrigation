% Define API Key and Location
api_key = 'aed4bcf416a982aca6a4d1d5df28cbde'; % Replace with your actual API key
location = 'Ghaziabad,IN'; % Change if needed
N = 48; % Forecast for 5 days (5 * 24 hours)

% Fetch Weather Forecast
[TA_forecast, RH_forecast, ~] = get_weather_forecast11(api_key, location, N);

% Calculate Dew Point Temperature
TD_forecast = TA_forecast - (100 - RH_forecast) / 5;

% Calculate Evaporation Rate
k = 0.12;
E_forecast = k * (TA_forecast - TD_forecast) .* (1 - RH_forecast / 100);

% Time Vector (in hours from current time)
time_hours = 0:N-1;

% Plot 1: Ambient Temperature and Humidity
figure;
yyaxis left
plot(time_hours, TA_forecast, 'r-', 'LineWidth', 2);
ylabel('Ambient Temperature (°C)');
hold on;
yyaxis right
plot(time_hours, RH_forecast, 'b-', 'LineWidth', 2);
ylabel('Relative Humidity (%)');
xlabel('Time (hours from now)');
title('Ambient Temperature and Humidity Forecast');
legend('Temperature', 'Humidity');
grid on;

% Plot 2: Evaporation Rate
figure;
plot(time_hours, E_forecast, 'g-', 'LineWidth', 2);
xlabel('Time (hours from now)');
ylabel('Evaporation Rate');
title('Evaporation Rate Forecast');
grid on;
