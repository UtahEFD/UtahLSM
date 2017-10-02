close all;clear all;clc;
initfiledir = pwd;

zo         = 0.15;                    %roughness height
zlevel     = 10;                      %Height of measurement (meters)
dt         = 0.0007;                  %Normalized time step
z_i        = zlevel;                  %Nondimensional height parameter (current one assumes 'square' domain)
TimeSteps  = 0;
T_scale    = 1;
Q_scale    = 1;
u_star      = 0.45;
densityAir = 1.204;
Cp_air     = 1005.0;

save([initfiledir,'/zo.ini'],'zo','-ascii','-tabs')
save([initfiledir,'/zlevel.ini'],'zlevel','-ascii','-tabs')

zGnd = [0.0;...
    0.005;...
    0.01;...
    0.02;...
    0.04;...
    0.06;...
    0.08;...
    0.12;...
    0.20;...
    0.30;...
    0.50;...
    1.00;...
    2.00]; %meters

% fid = fopen([initfiledir,'/zGnd.ini'],'w','l');
% fwrite(fid,zGnd(:),'double');
% fclose(fid);

save([initfiledir,'/zGnd.ini'],'zGnd','-ascii','-tabs')

numSoilLevels = length(zGnd); % 13

soilMatrix = zeros(numSoilLevels,1);

soilZarrayMeas = [0, 2, 4, 8, 12, 20, 30, 50]/100;

soilTempMeas = [292.2952,...
    292.3390,...
    292.9320,...
    292.7451,...
    292.8211,...
    292.0319,...
    291.1353,...
    289.2415];

soilTempMeas(end+1) = 283.15;
soilZarrayMeas(end+1) = 2;

soilTempInit = interp1(soilZarrayMeas,soilTempMeas,zGnd);

% initialize matrices of soil properties to zero
porosity     = soilMatrix;
satPotential = soilMatrix;  % saturated potential
satHydrCond  = soilMatrix;  % saturated hydraulic conductivity
soilExponent = soilMatrix;  % b in Clapp/Hornberger, non-dimensional
heatCapSoil  = soilMatrix;  % volumetric heat capacity of pure soil

% temperature and moisture profile
temperatureProfile = soilMatrix;
moistureProfile    = soilMatrix;

% albedo information for each x, y surface point
% albedoFlag == 1 => function form of albedo with moisture content and elevation angle of sun
% albedoFlag == 0 => constant albedo    
% albedoFlag pertains to entire domain
albedoFlag = 1;  
albedo(1)    = 0.33;
minAlbedo(1) = 0.13; % only used if albedoFlag=1

porosity(:) = [ 0.47,...
    0.47,...
    0.47,...
    0.47,...
    0.47,...
    0.47,...
    0.673,...
    0.863,...
    0.863,...
    0.863,...
    0.863,...
    0.863,...
    0.863 ]; % m^3/m^3

satPotential(:) = -[0.405,...
    0.405,...
    0.405,...
    0.405,...
    0.405,...
    0.405,...
    0.381,...
    0.356,...
    0.356,...
    0.356,...
    0.356,...
    0.356,...
    0.356]   / z_i; % meters

satHydrCond(:) = [1.3,...
    1.3,...
    1.3,...
    1.3,...
    1.3,...
    1.3,...
    4.65,...
    8.0,...
    8.0,...
    8.0,...
    8.0,...
    8.0,...
    8.0] *1e-6 / u_star;   % m/s

soilExponent(:) = [11.4,...
    11.4,...
    11.4,...
    11.4,...
    11.4,...
    11.4,...
    9.575,...
    7.75,...
    7.75,...
    7.75,...
    7.75,...
    7.75,...
    7.75];

heatCapSoil(:) = [1090,...
    1090,...
    1090,...
    1090,...
    1090,...
    1090,...
    965,...
    840,...
    840,...
    840,...
    840,...
    840,...
    840] *1000 /(densityAir*Cp_air);


temperatureProfile(:) = soilTempInit;
temperatureProfile(:) = [287.450,...
    287.250,...
    288.677,...
    289.777,...
    290.335,...
    291.340,...
    291.080,...
    291.181,...
    290.520,...
    289.958,...
    288.689,...
    285.350,...
    283.150];


moistureProfile(:) = [0.24,...
    0.247,...
    0.269,...
    0.274,...
    0.310,...
    0.330,...
    0.330,...
    0.360,...
    0.450,...
    0.470,...
    0.470,...
    0.570,...
    0.570] ;% + 0.01 *(1-zGnd'/2);
bb=5;

% fid = fopen([initfiledir,'/soilTypeParams.ini'],'w','l');
% fwrite(fid,porosity,'double');
% fwrite(fid,satPotential,'double');
% fwrite(fid,satHydrCond,'double');
% fwrite(fid,soilExponent,'double');
% fwrite(fid,heatCapSoil,'double');
% fclose(fid);

save([initfiledir,'/zGnd.ini'],'zGnd','-ascii','-tabs')
soilout=[porosity';satPotential';satHydrCond';soilExponent';heatCapSoil'];
save([initfiledir,'/soilTypeParams.ini'],'soilout','-ascii','-tabs')

% fid = fopen([initfiledir,'/soilTemperature.ini'],'w','l');
% fwrite(fid,temperatureProfile,'double');
% fclose(fid);
save([initfiledir,'/soilTemperature.ini'],'temperatureProfile','-ascii','-tabs')

% fid = fopen([initfiledir,'/soilMoisture.ini',],'w','l');
% fwrite(fid,moistureProfile,'double');
% fclose(fid);
save([initfiledir,'/soilMoisture.ini'],'moistureProfile','-ascii','-tabs')

% fid = fopen([initfiledir,'/albedo.ini'],'w','l');
% fwrite(fid,albedo,'double');fclose(fid);
save([initfiledir,'/albedo.ini'],'albedo','-ascii','-tabs')
if albedoFlag
    save([initfiledir,'/minAlbedo.ini'],'minAlbedo','-ascii','-tabs')
%     fid = fopen([initfiledir,'/minAlbedo.ini'],'w','l');
%     fwrite(fid,minAlbedo,'double');fclose(fid);
end
    
fprintf('%s','done!');