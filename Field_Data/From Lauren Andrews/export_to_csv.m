load('site_3_2012.mat')
load('FOXX_weather_2012.mat')
load('foxx_wx_2011.mat')
load('foxx_wx_2012.mat')
load('FOXX_weather_2011.mat')

stage_data = site_3_2012.stage_m;
stage_time = site_3_2012.UTC_DOY;

temperature = foxx_wx_2012.airtemp_C;
time = foxx_wx_2012.UTC_DOY;
incom_sw_watt_data = foxx_wx_2012.incom_sw_wattsperm2;


UTC_DOY = foxx_wx_2012.UTC_DOY;
record = foxx_wx_2012.record;
refl_sw_wattsperm2 = foxx_wx_2012.refl_sw_wattsperm2;
incom_sw_wattsperm2 = foxx_wx_2012.incom_sw_wattsperm2;
albedo = foxx_wx_2012.albedo;
windspeed_mpers = foxx_wx_2012.windspeed_mpers;
winddir_deg = foxx_wx_2012.winddir_deg;
airtemp_C = foxx_wx_2012.airtemp_C;
relhumidity = foxx_wx_2012.relhumidity;
dz_meters = foxx_wx_2012.dz_meters;



weather_foxx = array2table([UTC_DOY, record, refl_sw_wattsperm2, incom_sw_wattsperm2, albedo, windspeed_mpers, winddir_deg, airtemp_C, relhumidity, dz_meters],'VariableNames',{'UTC_DOY', 'record', 'refl_sw_wattsperm2', 'incom_sw_wattsperm2', 'albedo', 'windspeed_mpers', 'winddir_deg', 'airtemp_C', 'relhumidity', 'dz_meters'});
writetable(weather_foxx,'weather_foxx_2012.csv')

% figure(1)
% plot(stage_time,stage_data)
% 
% figure(2)
% plot(temperature_time,temperature_data)
% 
% figure(3)
% plot(temperature_time,incom_sw_watt_data)