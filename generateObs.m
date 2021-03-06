%   Generate Observatios For MRL model
%   The link size attributes are already computed for all pairs OD
%   --------
%   filename:   Name of file which stores the observations
%   x0:         Given parameters
%   ODpairs:    Matrix with all OD pairs
%   nbobsOD:    Number of generated obs each OD
%%
function ok = generateObs(filename, x0, ODpairs, nbobsOD)

    global incidenceFull; 
    global Mfull;
    global Ufull;
    global Op;
    global isLinkSizeInclusive;
    
    % Generate Obs
    % ----------------------------------------------------  
    mu = 1; % MU IS NORMALIZED TO ONE
    % Parameter for the radom term
    location = 0;
    scale = 1;
    euler = 0.577215665;
    rng('shuffle');
    [lastIndexNetworkState] = size(incidenceFull, 1);
    nbobs = size(ODpairs,1);
    % For the OD independence attributes
    IregularNetwork = incidenceFull(1:lastIndexNetworkState,1:lastIndexNetworkState);    
    dummy = lastIndexNetworkState + 1;
    % Obs = sparse(zeros(nbobs * nbobsOD, dummy));
    % Loop over all OD pairs
    for n = 1:nbobs
        n
        %% Generate paramters
        Op.x = [normrnd(x0(1),x0(2));x0(3:Op.n)];
        %----------------------------
        dest = ODpairs(n, 1);
        orig = ODpairs(n, 2);
        % get Link Size attributes
        % Get M and U
        Mfull = getM(Op.x, isLinkSizeInclusive); % matrix with exp utility for given beta
        M = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
        addColumn = Mfull(:,dest);
        M(:,lastIndexNetworkState+1) = addColumn;
        M(lastIndexNetworkState+1,:) = zeros(1,lastIndexNetworkState+1);
        [Z, expVokBool] = getExpV(M); % vector with value functions for given beta                                                                     
        if (expVokBool == 0)
            disp('The parameters are not fesible')
            ok = false;
            return; 
        end    
        V = log(Z);
        % Get Utility
        Ufull = getU(Op.x, isLinkSizeInclusive); % matrix of utility for given beta
        U = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
        addColumn = Ufull(:,dest);
        U(:,lastIndexNetworkState+1) = addColumn;
        U(lastIndexNetworkState+1,:) = zeros(1,lastIndexNetworkState+1);         
        addColumn = incidenceFull(:,dest);
        Incidence = IregularNetwork;
        Incidence(:,lastIndexNetworkState+1) = addColumn;
          
        % Now we have all utilities, time to generate the observations      
        for i = 1: nbobsOD
            Obs((n-1)*nbobsOD + i, 1) = dest;
            Obs((n-1)*nbobsOD + i, 2) = orig;
            k = orig;
            t = 3;
            while k ~= dummy
                ind = find(Incidence(k,:));
                nbInd = size(ind,2);
                bestUtilities = -1e6;
                for j = 1: nbInd
                    utility = U(k,ind(j)) + V(ind(j)) + random('ev',location,scale) - euler ;                   
                    if bestUtilities < utility
                        bestInd = ind(j);
                        bestUtilities = utility;
                    end
                end
                if bestInd ~= dummy
                    Obs((n-1)*nbobsOD + i, t) = bestInd;
                    t = t + 1;
                end
                k = bestInd;
            end
            Obs((n-1)*nbobsOD + i, t) = dest;
        end
    end
    % Write to file::
    [i,j,val] = find(Obs);
    data_dump = [i,j,val];
    save(filename,'data_dump','-ascii');
    ok = true;
end