function [TA_forecast, RH_forecast, rain_forecast] = get_weather_forecast(api_key, location, N) 
    if nargin < 2 || isempty(location)
        location = 'Ghaziabad,IN';
        api_key = 'd6bc4ca03f710c49bb0b20df9e686187';
    end
    if nargin < 3
        N = 24;
    end

    url = sprintf('https://api.openweathermap.org/data/2.5/forecast?q=%s&units=metric&appid=%s',location, api_key);
    try
        data = webread(url);
    catch
        error('Weather API request failed. Check internet connection and API key.');
    end

    num_forecasts = min(N, length(data.list));
    TA_forecast = zeros(1, num_forecasts);
    RH_forecast = zeros(1, num_forecasts);
    rain_forecast = zeros(1, num_forecasts);

    for i = 1:num_forecasts
        TA_forecast(i) = data.list(i).main.temp;
        RH_forecast(i) = data.list(i).main.humidity;
        if isfield(data.list(i), 'rain') && isfield(data.list(i).rain, '3h')
            rain = data.list(i).rain.('3h') / 3;
        else
            rain = 0;
        end
        rain_forecast(i) = rain;
    end

    forecast_times = 0:3:3*(num_forecasts-1);
    query_times = 0:N-1;

    if length(forecast_times) < 2
        error('Insufficient data for interpolation.');
    end

    TA_forecast = interp1(forecast_times, TA_forecast, query_times, 'spline', 'extrap');
    RH_forecast = interp1(forecast_times, RH_forecast, query_times, 'spline', 'extrap');
    rain_forecast = interp1(forecast_times, rain_forecast, query_times, 'spline', 'extrap');

    TA_forecast = max(TA_forecast, 0);
    RH_forecast = max(min(RH_forecast, 100), 0);
    rain_forecast = max(rain_forecast, 0);
end
