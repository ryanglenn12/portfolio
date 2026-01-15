% original DAQ generator below - now use Transect script
% %% DATA STRUCTURE CREATION
% % I added this section in to create data structure of CO concs with time
% % used later to calculate dilution rate constant
% time = [0; 3240; 3420; 3720; 4020; 4620; 4740; 6960; 7140];
% CO = [132; 195; 149; 168; 175; 205; 250; 220; 245];
% dataStruc.time = time;
% dataStruc.CO = CO;
% save('C:\Users\sarp\Projects\SARP\F0AM-4.3.0.1\Docs', 'dataStruc');

%% TRANSECT DATA GENERATION (Standalone)
% call Palisades_data file 
Palisades_data

%% METEOROLOGY
%{
Pressure, temperature, and either RH or H2O are required Met inputs.
Dilution is calculated using the change in CO over each model step.
All calculated J-values will be scaled to J4.
%}

%kdil calculation using CO decay and inversion of dilution equation
% dX/dt = -kdil*(X - Xb)
dCOdt        = diff(DAQ.CO)./diff(DAQ.TIME); %loss rate
dCOdt(end+1) = dCOdt(end);
COmid        = (DAQ.CO + [DAQ.CO(2:end);DAQ.CO(end)])/2; %CO in middle of step
kdil         = -dCOdt./(COmid-95);

% similar calculation for tgauss
%dX/dt = (-1/(tgauss + 2t))*(X - Xb)
tgauss = 1./kdil - 2.*DAQ.TIME;
% % doesn't work great (gives negative numbers, probably b/c not in center of plume)
% % instead just fit it to kdil for best guess
% tgauss = fminsearch(@(x) sum((kdil - 1./(x + 2*DAQ.TIME)).^2),300);

% USER: try running with either tgauss or kdil and see how results compare!

Met = {...
    'P'             DAQ.P;... %Pressure, mbar
    'T'             DAQ.T;...    %Temperature, K
    'RH'            DAQ.RH;...    %Relative Humidity, percent    
    % 'kdil'          kdil;...    %Dilution constant, /sec     
    'tgauss'        tgauss; %alternative Gaussian dilution initial timescale
    'SZA'           DAQ.SZA;... %solar zenith angle
    'J4'            DAQ.JNO2;... %NO2 photolysis frequency;
    'jcorr'         {'J4'};... %correction factor for NO2 and O3 because of its capacity 
    % to oxidize and make ozone, plus literature supporting how often it is
    'ALT'           600;...  % altitude, meters
    'O3col'         300;...  % ozone column 
    'albedo'        0.075;... %albedo between 0-1; 0.1 for typical vegetation
    % 'LFlux'     'PalisadesLightFlux.txt'
    };

clear dCOdt COmid kdil tgauss

%% CHEMICAL CONCENTRATIONS
%{
The first observational point is used as initial inputs. 
All concentrations will be calculated "free running," meaning no constraints along the transect.
%}
InitConc = {...
%   names               conc(ppb)               HoldMe    
    'O3'                DAQ.O3(1)               0;... 
    'NO2'               DAQ.NO2(1)              0;... 
    'CH4'               DAQ.CH4(1)              0;...
    'CO'                DAQ.CO(1)               0;...
    % 'CO2'               DAQ.CO2(1)              0;...
    'NO'                DAQ.NO(1)               0;...

    'CH3OH'             DAQ.methanol(1)         0;... %Methanol
    'C2H5OH'            DAQ.ethanol(1)         0;...  %ethanol
    'CH3CHO'            DAQ.acetaldehyde(1)    0;... %acetaldehyde
    'C3H6'              DAQ.propene(1)          0;... %propene
    'BENZENE'           DAQ.benzene(1)          0;... %benzene
    'FURAN'             DAQ.furan(1)            0;... %furan; only decay!
    'C5H8'              DAQ.isoprene(1)         0;... %isoprene
    'CH3COCH3'          DAQ.acetone(1)          0;... %acetone
    'C2H2'              DAQ.ethyne(1)           0;... %ethyne
    'C2H6'              DAQ.ethane(1)           0;... %Ethane
    'C2H4'              DAQ.ethene(1)           0;... %Ethene
    'C3H8'              DAQ.propane(1)          0;... %propane
    'NC4H10'            DAQ.nbutane(1)          0;... %n-butane
    'IC4H10'            DAQ.ibutane(1)          0;... %i-butane
    'IC5H12'            DAQ.ipentane(1)         0;... %i-pentane
    'NC5H12'            DAQ.npentane(1)         0;... %n-pentane
    'TOLUENE'           DAQ.toluene(1)         0;... %toluene
    'CH3CL'             DAQ.dichloromethane(1) 0;... % chloromethane
    'CH2CL2'            DAQ.chloromethane(1)   0;... % dichloromethane
    };

VOC_species = {'CH3OH', 'C2H5OH', 'CH3CHO', 'C3H6', 'BENZENE', 'FURAN',...
    'C5H8', 'CH3COCH3', 'C2H2', 'C2H6', 'C2H4', 'C3H8', 'NC4H10', 'IC4H10', ...
    'IC5H12', 'NC5H12', 'TOLUENE', 'CH3CL', 'CH2CL2'};

%% CHEMISTRY
%{
The ChemFiles input is a cell array of strings specifying functions and scripts for the chemical mechanism.
THE FIRST CELL is always a function for generic K-values.
THE SECOND CELL is always a function for J-values (photolysis frequencies).
All other inputs are scripts for mechanisms and sub-mechanisms.
Here we give an example using MCMv3.3.1.  Note that this mechanism was extracted from the MCM website for
the specific set of initial species included above.
"FURFURAL_FURAN" is a very simple set of reactions for initial oxidation of these species,
which are not included in MCM. For extra fun, try toggling this on and off to compare results.
%}

ChemFiles = {...
    'MCMv331_K(Met)';
    'MCMv331_J(Met,2)'; %Jmethod flag of 0 specifies "MCM" J-value method, 1 bottom up, J hybrid
    'MCMv331_LGPlumeSubset';
    'FURFURAL_FURAN'; %very simple initial oxidation
    };
 
%% DILUTION CONCENTRATIONS
% Background concentrations are taken from observations just outside the plume.

BkgdConc = {...
%   names               values   %0 for all zeros, 1 to use InitConc
    'O3'                background_O3               ;... 
    'NO2'               background_NO2              ;... 
    'CH4'               background_CH4              ;...
    'CO'                background_CO               ;...
    % 'CO2'               background_CO2              ;...
    'NO'                background_NO               ;...

    'CH3OH'             background_methanol         ;... %Methanol
    'C2H5OH'            background_ethanol          ;...  %ethanol
    'CH3CHO'            background_acetaldehyde     ;... %acetaldehyde
    'C3H6'              background_propene          ;... %propene
    'BENZENE'           background_benzene          ;... %benzene
    'FURAN'             background_furan            ;... %furan; only decay!
    'C5H8'              background_isoprene         ;... %isoprene
    'CH3COCH3'          background_acetone          ;... %acetone
    'C2H2'              background_ethyne           ;... %ethyne
    'C2H6'              background_ethane           ;... %Ethane
    'C2H4'              background_ethene           ;... %Ethene
    'C3H8'              background_propane          ;... %propane
    'NC4H10'            background_nbutane          ;... %n-butane
    'IC4H10'            background_ibutane          ;... %i-butane
    'IC5H12'            background_ipentane         ;... %i-pentane
    'NC5H12'            background_npentane         ;... %n-pentane
    'TOLUENE'           background_toluene          ;... %toluene
    'CH3CL'             background_dichloromethane  ;... % chloromethane
    'CH2CL2'            background_chloromethane    ;... % dichloromethane
    };

%% %% OPTIONS
%{
"Verbose" can be set from 0-3; this just affects the level of detail printed to the command
  window regarding model progress.
"EndPointsOnly" is set to 0 because we want output to include all concentrations along each model step.
"LinkSteps" is set to 1 because each step is connected.
"IntTime" is the integration time for each step. The average for constraints is 250s.
"SavePath" is just a filename, which will be saved by default in the \Runs\ folder.
%}

ModelOptions.Verbose = 1;         %flag for verbose output
ModelOptions.EndPointsOnly = 0;   %flag for concentration and rate outputs
ModelOptions.LinkSteps = 1;       %flag for using end-points of one run to initialize next run
ModelOptions.IntTime = 3600;       %integration time for each step, seconds
ModelOptions.SavePath = 'Palisades804Results';
ModelOptions.GoParallel    = 0;

%% MODEL RUN
% Now we call the model.
% Output will be saved in the "SavePath" above and will also be written to the structure S.
% Let's also throw away the inputs (don't worry, they are saved in the output structure).

S = F0AM_ModelCore(Met,InitConc,ChemFiles,BkgdConc,ModelOptions);
clear Met InitConc ChemFiles BkgdConc ModelOptions

%% FIGURES AND ANALYSIS

% calculate normalized excess mixing ratios: delta_X_dil = delta_X * delta_CO_source/delta_CO
% this is standard procedure for biomass burning plumes.
% We will also put these into a new structure compatible with F0AM plotting routines.
EMR_CO = S.Conc.CO - S.BkgdConc.CO(1);
NEMR_CO = EMR_CO(1)./EMR_CO;

Sd.Met = S.Met; Sd.Cnames = S.Cnames; Sd.Time = S.Time; Sd.BkgdConc = S.BkgdConc;
for i=1:length(S.Cnames)
    name = S.Cnames{i};
    if isfield(S.BkgdConc,name), b = S.BkgdConc.(name)(1);
    else b = 0;
    end
    Sd.Conc.(name) = (S.Conc.(name) - b).*NEMR_CO;
end

% % First, check that dilution is working as it should
% PlotConc('CO',S);
% hold on
% plot(DAQ.TIME, DAQ.CO,'k+','markersize',18,'linewidth',4)
% 
% % Now let's look at NOy speciation
% % first we have to get all NOy, which is a tall order for MCM, but we can try by using 
% % SMILES strings to identify functional groups.
% 
% PNs = SearchSMILES('peroxyNitrate',S.Cnames,'v331'); %peroxy nitrates
% ANs = SearchSMILES('alkylNitrate',S.Cnames,'v331'); %alkyl nitrates
% otherNOy = {'NO';'NO2';'NO3';'N2O5';'HONO';'HNO3';'HO2NO2'};
% [~,iother] = ismember(otherNOy,S.Cnames);
% iNOy = unique([iother; PNs.index; ANs.index]);
% 
% NOygroups = {{'NOx','NO','NO2'},['PNs';S.Cnames(PNs.index)],['ANs';S.Cnames(ANs.index)],'HONO','HNO3'};
% PlotConcGroup(NOygroups,Sd,5,'sortem',0,'name','NO_y')
% 
% %Next, lets look at ozone production.
% %modeled vs observed
% 
% EMR_O3 = Sd.Conc.O3 - Sd.BkgdConc.O3(1); % if you have background O3
% NEMR_O3 = EMR_O3/EMR_CO;
% figure
% PlotConc('O3', S);
% hold on
% plot(DAQ.TIME, DAQ.O3,'g','markersize',18,'linewidth',4)
% 
% figure
% plot(Sd.Time, Sd.Conc.O3 - Sd.BkgdConc.O3(1), 'b-',...
%     Sd.Time, EMR_O3, 'm-')
% xlabel('Reaction Time (s)')
% ylabel('Ozone Mixing Ratios (ppb)')
% legend('Modeled', 'Observed')
% 
% figure
% plot(Sd.Time, EMR_O3, 'ro')
% xlabel('Reaction Time (s)')
% ylabel('Modeled Ozone Production (ppb)')

figure
plot(Sd.Time, Sd.Conc.O3, 'b-')
xlabel('Reaction Time (s)')
ylabel('Normalized Excess Mixing Ratio of Ozone')


% % It is easiest to look at production of ozone's precursor, NO2.
% % Also, since dilution is so strong, we will look at this in a rate-normalized sense
% % because we are just interested in speciation.
% PlotRates('NO2',S,5,'sumEq',1,'scale',0)
% 
% PlotRates('NO',S,5,'sumEq',1,'scale',0)
% 
% % it looks like H02 is the dominant source after NO photolysis. Let's keep going.
% % Note that FURFURAL, which is not in the MCM base chemistry, is a non-trivial source of HO2.
% PlotRates('HO2',S,5,'sumEq',1,'scale',0)
% 
% % DILUTION IN NEGATIVES AT BEGINNING FOR SOME REASON?
% PlotRates('HCHO',S,5,'sumEq',1,'scale',0)
% 
% % encountering errors for following graphs bc axes labels and graph styles
% % are tranferring over, so specify for each graph, and find a way to clear
% % previous graph settings
% %same thing as above but assuming not much change is happening with time
% PlotRatesAvg('HO2', S, 5)
% PlotRatesAvg('HCHO', S, 5)
% 
% 
% PlotRatesAvg('CH3COCH3', S, 5)
% PlotRatesAvg('CH3CHO', S, 5)
% % plots two indicator ratios for NOX and VOCs based on literature
% ratio1 = Sd.Conc.H2O2./Sd.Conc.HNO3;
% ratio2 = Sd.Conc.HCHO./Sd.Conc.NO2;
% figure
% plot(Sd.Time, ratio1, 'g')
% xlabel('Reaction Time (s)')
% ylabel('H2O2:HNO3')
% 
% figure
% plot(Sd.Time, ratio2, 'g')
% xlabel('Reaction Time (s)')
% ylabel('HCHO:NO2')
% 
% % It looks like speciation doesn't change much over the integration period, so we can continue
% % investigation in an average sense. 
% % CO and HCHO are signficant sources of HO2. Most of the CO is primary emissions, and so is HCHO
% % (intitally). What about secondary HCHO?
% PlotRatesAvg('HCHO',S,5)
% PlotRatesAvg('CH3O2',S,5,'sumEq',1)
% PlotRatesAvg('CH3CO3',S,5,'sumEq',1)
% % Secondary formaldehyde is mostly from acetaldehyde. Note this will changes as the plume ages more.
