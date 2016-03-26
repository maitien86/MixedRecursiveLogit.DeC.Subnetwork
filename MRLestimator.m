%   MAIN
%% Parameter: 1. number of draws, 2. real or sml 3. LS or noLS 4. IRN or CRN 5. RP or EC
function reTxt = MRLestimator(ndraws, isLS, mixedtype, corrtype, numberOfEC)    
    Credits;
    global Gradient;
    global Op;
    global file_AttEstimatedtime;
    global file_linkIncidence;
    global file_observations;
    global file_turnAngles;
    global isLinkSizeInclusive;
    global nbobs;  
    global isFixedUturn;
    global nDraws;
    global mixedType;
    global corrType;
    global resultsTXT; 
    global nECs;

    
    %% Set number of draws
    if isnumeric(ndraws)
        nDraws = ndraws;
    else
        nDraws = str2double(ndraws);
    end    
    %% Set number of Error components
    if isnumeric(numberOfEC)
        nECs = numberOfEC;
    else
        nECs = str2double(numberOfEC);
    end
    
    %% Data files
    file_linkIncidence = './Input/linkIncidence.txt';
    file_AttEstimatedtime = './Input/TravelTime.txt';
    file_turnAngles = './Input/TurnAngles.txt';
    RealObs = './simulatedData/ObsGLS.txt';
    
    % isLS: with global LS attribute
    % noLS: without global LS attribute
    if strcmp(isLS,'LS')
        isLinkSizeInclusive = true;
    else
        isLinkSizeInclusive = false;
    end    
    file_observations = RealObs;        
    
    global Omega;
    rng('shuffle');

    % mixedtype = 'IRN' or 'CRN'
    if strcmp(mixedtype,'IRN')
        mixedType = OptimizeConstant.IRN;
        Omega = randn(nbobs, nDraws);
    else
        mixedType = OptimizeConstant.CRN;
        if strcmp(corrtype,'RP')
            corrType = OptimizeConstant.RP;
            Omega = randn(nDraws, 1);
        else
            corrType = OptimizeConstant.EC;
            Omega = randn(nDraws, nECs);
        end
    end
    
    %% Initialize the optimizer structure
    isFixedUturn = false;
    loadData;
    %% Load subnetworks
    if corrType == OptimizeConstant.EC
        LoadSubNetworks;
    end

    %% Create optimization structure
    Op = Op_structure;
    initialize_optimization_structure();
    Op.x(Op.natt + 1: Op. n) = Op.x(Op.natt + 1: Op. n) * 0;
    %% Load attribute to matrices
    getAtt();
    Gradient = zeros(nbobs,Op.n);
    % For EC 
    if nECs == 5
        if isLinkSizeInclusive == false
            %Op.x = [-2.715312e+00; -9.270908e-01;-4.089189e-01;-4.428802e+00;-1.138836e+00;-5.479881e-03;5.661129e-01;1.715967e+00;4.793310e-01];
            Op.x = [-2.715312e+00; -9.270908e-01;-4.089189e-01;-4.428802e+00;0;0;0;0;0];
        else
            %Op.x = [-2.763899e+00; -9.330285e-01;-3.675629e-01;-4.455489e+00;-3.205717e+00;-1.267890e+00;5.277898e-02;1.600985e-01;-1.507828e+00;-1.077993e+00];
            Op.x = [-2.763899e+00; -9.330285e-01;-3.675629e-01;-4.455489e+00;-3.205717e+00;0;0;0;0;0];
        end
    else
        if isLinkSizeInclusive == false
            %Op.x = [-2.715312e+00; -9.270908e-01;-4.089189e-01;-4.428802e+00;-1.138836e+00;-5.479881e-03;5.661129e-01;1.715967e+00;4.793310e-01];
            Op.x = [-2.715312e+00; -9.270908e-01;-4.089189e-01;-4.428802e+00;0;0;0;0;0;0];
        else
            %Op.x = [-2.763899e+00; -9.330285e-01;-3.675629e-01;-4.455489e+00;-3.205717e+00;-1.267890e+00;5.277898e-02;1.600985e-01;-1.507828e+00;-1.077993e+00];
            Op.x = [-2.763899e+00; -9.330285e-01;-3.675629e-01;-4.455489e+00;-3.205717e+00;0;0;0;0;0;0];
        end
    end
    
    
    initialize_switching_structure();
    Op.Switching = OptimizeConstant.SW_RETRO;
    Op.Optim_Method = OptimizeConstant.TRUST_REGION_METHOD;
    Op.Hessian_approx = OptimizeConstant.BFGS;
    %% Starting optimization
    tic ;
    disp('Start Optimizing ....')
    [Op.value, Op.grad ] = LL(Op.x);
    Op.delta = 0.1 * norm(Op.grad);
    PrintOut(Op);
    %% ResultsTXT

    header = [sprintf('%s \n',file_observations) Op.Optim_Method];
    header = [header sprintf('\nNumber of observations = %d \n', nbobs)];
    header = [header sprintf('Hessian approx methods = %s \n', OptimizeConstant.getHessianApprox(Op.Hessian_approx))];
    header = [header sprintf('Mixed method = %s - %s \n',  corrtype, numberOfEC)];   
    resultsTXT = header;
    
    %% LOOP 
    while (true)    
      Op.k = Op.k + 1;
      if strcmp(Op.Optim_Method,OptimizeConstant.LINE_SEARCH_METHOD);
        ok = line_search_iterate();
        if ok == true
            PrintOut(Op);
        else
            disp(' Unsuccessful process ...')
            break;
        end
      else
        btr_interate();
        %btr_swretro_interate();
        PrintOut(Op);
      end
      [isStop, Stoppingtype, isSuccess] = CheckStopping(Op);  
      % Check stopping criteria
      if(isStop == true)
          isSuccess
          fprintf('The algorithm stops, due to %s', Stoppingtype);
          resultsTXT = [resultsTXT sprintf('The algorithm stops, due to %s \n', Stoppingtype)];
          break;
      end
      
    end
    %%   Compute variance - Covariance matrix
 %%   Compute variance - Covariance matrix
    resultsTXT = header;
    PrintOut(Op);
    fprintf('\n Calculating VAR-COV ...\n');
   
    %%
    getCov;
    if mixedType == OptimizeConstant.IRN
        method = 'IRN';
    else 
        method = 'CRN';
    end
    
    %% Finishing, send email notification ...
    ElapsedTtime = toc
    resultsTXT = [resultsTXT sprintf('\n Number of function evaluation %d \n', Op.nFev)];
    resultsTXT = [resultsTXT sprintf('\n Estimated time %d \n ----------------------------------------------------------------- \n', ElapsedTtime)];
    
    fileID = fopen('MixedResults_realdataEC.txt','at+');
    fprintf(fileID,resultsTXT);
    fclose(fileID);
end
