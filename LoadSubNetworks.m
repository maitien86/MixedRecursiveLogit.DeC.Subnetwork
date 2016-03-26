%% Load subnetworks
sub = objArray(5);
sub(1).value = spconvert(load('./Subnetworks/R50N_mapped.txt'))';
sub(2).value = spconvert(load('./Subnetworks/R50S_mapped.txt'))';
sub(3).value = spconvert(load('./Subnetworks/R70N_mapped.txt'))';
sub(4).value = spconvert(load('./Subnetworks/R70S_mapped.txt'))';
sub(5).value = spconvert(load('./Subnetworks/RC_mapped.txt'))';

%% Mapping travel times
global EstimatedTime;
global TTSub;
TTSub = objArray(6);
maxState = size(EstimatedTime,2);
B = 1:maxState;
for i = 1:5
    A = ismember(B,sub(i).value);
    TTSub(i).value = (bsxfun(@times,A,EstimatedTime)); 
end
TTSub(6).value = EstimatedTime;
