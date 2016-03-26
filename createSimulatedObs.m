%   Generate observations for MRL model
%%
globalVar;

file_linkIncidence = './Input/linkIncidence.txt';
file_AttEstimatedtime = './Input/ATTRIBUTEestimatedtime.txt';
file_turnAngles = './Input/ATTRIBUTEturnangles.txt';
file_observations = './Input/observationsForEstimBAI.txt';

% without LS
% 
isLinkSizeInclusive = false;
isFixedUturn = false;
loadData;
Op = Op_structure;
initialize_optimization_structure();
%Op.n = 4;
disp('Observation is generating ....')
x0 = [-2.0, 0.6, -1.0,-1.0,-4.0]';
ODpairs = Obs(:,1:2);
nbobsOD = 1;
filename = './simulatedData/Obs_NoGLS.0.6.txt';
generateObs(filename, x0, ODpairs, nbobsOD);


% with LS

isLinkSizeInclusive = true;
isFixedUturn = false;
loadData;
Op = Op_structure;
initialize_optimization_structure();
%Op.n = 4;
disp('Observation is generating ....')
x0 = [-2.0, 0.6, -1.0,-1.0,-4.0,-3.0]';
ODpairs = Obs(:,1:2);
nbobsOD = 1;
filename = './simulatedData/Obs_GLS.0.6.txt';
generateObs(filename, x0, ODpairs, nbobsOD);

