% Palisades data generation script

% Define input sets for Gaussian dilution model
n_points = 30; % Number of points along the transect (e.g., 0 to 10 km in 0.5 km steps)
distance_km = linspace(1, 30, n_points)'; % Distance from source (km), transposed to column
wind_speed = 25; % Wind speed (m/s)
time_seconds = distance_km * 1000 / wind_speed; % Lagrangian time (s), assuming linear advection

% Species initial and background concentrations (ppb)
% bkgd concs taken from min values of California Air Resources Board 2021
% report of RECAP-CA airborne campaign
% init concs taken from avg values of late morning section of sampling data
% can alternatively break sampling into 2 groups by time: one bkgd, one init
% GHGs and 15 most concentrated VOCs from sampling data plus furan
initial_CO = 219;  background_CO = 119;
% initial_CO2 = 444; background_CO2 = 414875;
initial_CH4 = 2115;  background_CH4 = 2033;
initial_O3 = 34;  background_O3 = 40;           % is 40 too high a bkgd conc?
initial_NO2 = 37.6;  background_NO2 = 8.5;
initial_NO = 31.8; background_NO = 1.5;

initial_ethane = 3.48;  background_ethane = 2.27;  
initial_methanol = 2.86;  background_methanol = 2.17; 
initial_acetone = 3.20;  background_acetone = 2.55;
initial_propane = 1.54;  background_propane = 0.93;
initial_ethanol = 1.27;  background_ethanol = 0.67;
initial_acetaldehyde = 1.10;  background_acetaldehyde = 0.93;
initial_nbutane = 0.72;  background_nbutane = 0.36;
initial_ethene = 0.60;  background_ethene = 0.293;
initial_ethyne = 0.45;  background_ethyne = 0.272;
initial_ipentane = 0.29; background_ipentane = 0.136;
initial_ibutane = 0.29;  background_ibutane = 0.14;
initial_toluene = 0.16;  background_toluene = 0.07;
initial_propene = 0.15;  background_propene = 0.068;
initial_benzene = 0.13;  background_benzene = 0.074;
initial_npentane = 0.13;  background_npentane = 0.066;
initial_isoprene = 0.05;  background_isoprene = 0.16;
initial_furan = 0.018;  background_furan = 0.0035;
initial_dichloromethane = 0.120;  background_dichloromethane = 0.08;
initial_chloromethane = 0.574;  background_chloromethane = 0.545;

% Gaussian dilution model
tgauss = 2000; % was at 107 
CO_conc = initial_CO * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_CO;
% CO2_conc = initial_CO2 * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_CO2;
CH4_conc = initial_CH4 * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_CH4;
O3_conc = initial_O3 * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_O3;
NO2_conc = initial_NO2 * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_NO2;
NO_conc = initial_NO * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_NO;

ethane_conc = initial_ethane * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_ethane;
methanol_conc = initial_methanol * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_methanol;
acetone_conc = initial_acetone * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_acetone;
propane_conc = initial_propane * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_propane;
ethanol_conc = initial_ethanol * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_ethanol;
acetaldehyde_conc = initial_acetaldehyde * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_acetaldehyde;
nbutane_conc = initial_nbutane * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_nbutane;
ethene_conc = initial_ethene * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_ethene;
ethyne_conc = initial_ethyne * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_ethyne;
ipentane_conc = initial_ipentane * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_ipentane;
ibutane_conc = initial_ibutane * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_ibutane;
toluene_conc = initial_toluene * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_toluene;
propene_conc = initial_propene * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_propene;
benzene_conc = initial_benzene * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_benzene;
npentane_conc = initial_npentane * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_npentane;
isoprene_conc = initial_isoprene * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_isoprene;
furan_conc = initial_furan * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_furan;
dichloromethane_conc = initial_dichloromethane * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_dichloromethane;
chloromethane_conc = initial_chloromethane * exp(-time_seconds ./ (tgauss + 2 * time_seconds)) + background_chloromethane;

% Synthetic meteorological data
pressure_mbar = 880 * ones(n_points, 1); % Constant pressure (mbar), 21x1
temperature_K = 288 * ones(n_points, 1);  % Constant temperature (K), 21x1
relative_humidity = 10 * ones(n_points, 1); % Constant RH (%), 21x1
solar_zenith_angle = 55.26 * ones(n_points, 1); % Constant SZA (degrees), 21x1
JNO2 = 0.006 * ones(n_points, 1); % Constant NO2 photolysis frequency (1/s), 21x1

% Create DAQ structure with individual 21x1 numeric arrays
DAQ.TIME = time_seconds;
DAQ.CO = CO_conc;
% DAQ.CO2 = CO2_conc;
DAQ.CH4 = CH4_conc;
DAQ.O3 = O3_conc;
DAQ.NO2 = NO2_conc;
DAQ.NO = NO_conc;
DAQ.P = pressure_mbar;
DAQ.T = temperature_K;
DAQ.RH = relative_humidity;
DAQ.SZA = solar_zenith_angle;
DAQ.JNO2 = JNO2;

DAQ.ethane = ethane_conc;
DAQ.methanol = methanol_conc;
DAQ.acetone = acetone_conc;
DAQ.propane = propane_conc;
DAQ.ethanol = ethanol_conc;
DAQ.acetaldehyde = acetaldehyde_conc;
DAQ.nbutane = nbutane_conc;
DAQ.ethene = ethene_conc;
DAQ.ethyne = ethyne_conc;
DAQ.ipentane = ipentane_conc;
DAQ.ibutane = ibutane_conc;
DAQ.toluene = toluene_conc; 
DAQ.propene = propene_conc; 
DAQ.benzene = benzene_conc;
DAQ.npentane = npentane_conc;
DAQ.isoprene = isoprene_conc;
DAQ.furan = furan_conc;
DAQ.dichloromethane = dichloromethane_conc;
DAQ.chloromethane = chloromethane_conc;

% % plot ozone concentrations with time
% figure
% plot(DAQ.TIME, O3_conc)

% Save DAQ structure to out_file.m
output_file = fullfile('C:\Users\F0AM-4.3.0.1\Docs', 'Palisades_DAQ1.m');
fid = fopen(output_file, 'w');
if fid == -1
    error('Could not open %s for writing.', output_file);
end
fprintf(fid, '%% Generated DAQ structure from ExampleSetup_LagrangianPlume.m\n');
fprintf(fid, '%% Date: %s\n', '10:30 PM PDT on Thursday, July 24, 2025');
fprintf(fid, 'DAQ = struct();\n');
fields = fieldnames(DAQ);
for i = 1:length(fields)
    field = fields{i};
    data = DAQ.(field);
    
    if isnumeric(data) && isvector(data)
        fprintf(fid, 'DAQ.%s = [%.6g', field);
        fprintf(fid, ' %.6g', data(2:end));
        fprintf(fid, '];\n');
    else
        fprintf('⚠️ Field "%s" is not numeric — skipping\n', field);
    end
end

fclose(fid);
fprintf('DAQ structure saved as individual 21x1 numeric arrays in %s\n', output_file);
