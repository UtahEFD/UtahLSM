close all;clear all;clc;
initFileDir = pwd;

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

numSoilLevels = length(zGnd);

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
    0.356]; % meters

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
    8.0] *1e-6;   % m/s

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
    840] *1000;


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
    0.570];

soilout=[porosity';satPotential';satHydrCond';soilExponent';heatCapSoil'];

% write input files
save([initFileDir,'/soilLevels.ini'],'zGnd','-ascii','-tabs')
save([initFileDir,'/soilTypeParams.ini'],'soilout','-ascii','-tabs')
save([initFileDir,'/soilTemperature.ini'],'temperatureProfile','-ascii','-tabs')
save([initFileDir,'/soilMoisture.ini'],'moistureProfile','-ascii','-tabs')
    
fprintf('%s','done!');