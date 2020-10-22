load surfacedischarge.mat
%%
UTC_DOY = surface_discharge.UTC_DOY';
UTC_DOY_1h = surface_discharge.UTC_DOY_1h';
m3_m3s_6h = surface_discharge.m3_m3s_6h';
m4_m3s_6h = surface_discharge.m4_m3s_6h';
mF_m3s_6h = surface_discharge.mF_m3s_6h';
m3_m3s_1h_24hS = surface_discharge.m3_m3s_1h_24hS;
m4_m3s_1h_24hS = surface_discharge.m4_m3s_1h_24hS;
mF_m3s_1h_24hS = surface_discharge.mF_m3s_1h_24hS;
foxx_z_meters = surface_discharge.foxx_z_meters';
foxx_dz_mper6h = surface_discharge.foxx_dz_mper6h';

UTC_SOY = UTC_DOY*3600*24;
UTC_SOY_1h = UTC_DOY_1h*3600*24;

table6h = array2table([UTC_DOY,UTC_SOY,foxx_z_meters,foxx_dz_mper6h,m3_m3s_6h,m4_m3s_6h,mF_m3s_6h],'VariableNames', {'UTC_DOY','UTC_SOY','foxx_z_meters','foxx_dz_mper6h','m3_m3s_6h','m4_m3s_6h','mF_m3s_6h'});
table1h = array2table([UTC_DOY_1h,UTC_SOY_1h,m3_m3s_1h_24hS,m4_m3s_1h_24hS,mF_m3s_1h_24hS],'VariableNames', {'UTC_DOY','UTC_SOY_1h', 'm3_m3s_1h_24hS','m4_m3s_1h_24hS','mF_m3s_1h_24hS'});

writetable(table6h,'surface_discharge_andrews2014_6h.csv')
writetable(table1h,'surface_discharge_andrews2014_1h_smooth.csv')