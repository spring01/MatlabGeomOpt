
% Training DAgger with the expert with step size of 2
% Input:  "Pure Linear Combination" Policy Parameters (the mat file generated from PureLinearCombination.m)
% Output: Policy parameters (Coefficient: allCs)

clear classes;
allFits = AllFits;

%
nIteration = 25;

currentInd = 0;
step_size=2;
outputCell={};

%outputCell{1, HFiter} store the expert's demonstraion [dn  (ps)]
%targetEnergy          store the expert's demonstraion E

% The related parameters are stored in the mat file
% Including: Path, FragRange, FragRangeTest, envs
load('PureLinearCombinationResult.mat')
warning('off','all')

oriCs=allCs;
allCs={};

% The very initial
DCM = DistributeCoeff(Path,FragRange,envs);
%DCM.errorType=2;
DCM.errorType=0;

% from d0 -> d1
[err,Jac] = DCM.errJac(1);
% this is d1
Ps_iter = DCM.currentPs;
outputCell{1, 1} = Ps_iter;

% Refer to DAgger algorithm
expertDecay = 0.5;

for j=2:HF_Iteration/step_size
    %outputCell{1, HFiter} set the expert's demonstraion [dn  (ps)]
    outputCell{1, j} = recordP{(j-1)*step_size+1};
    disp([num2str(j) ' ==>' num2str((j-1)*step_size+1) ]);
end

for HFiter = 2:HF_Iteration/step_size 

    HFiter
    disp(['== HFiter: ' num2str(HFiter) ' ==']);
    disp([num2str((HFiter-2)*step_size+1) ' ==>' num2str((HFiter-1)*step_size+1) ]);

    % We set all zero as our initial guess in the beginning of each HFiter
    %lastX = zeros(size(ltrain.pars));
    
    %Cs = zeros(HFiter,1);
    %Cs(end-1:end) = [0 1]
    Cs = oriCs{HFiter}
    
    startInd = currentInd +1;
    
    DCM.objectiveDensities = recordP{(HFiter-1)*step_size+1};
    DCM.objectiveEs =        recordE{(HFiter-1)*step_size+1};

    for DaggerIter=1:HFiter-1
        
        currentInputStr=num2str(DaggerIter);
        expertRatio = expertDecay^(DaggerIter-1)
        
        
        ExpertDemonstration=cell(size(DCM.Ps));
        ExpertDemonstration(1,:,:)=DCM.Ps(1,:,:);
        for i=2:HFiter
            ExpertDemonstration(i,:,:) = outputCell{1,i-1} ;
        end

        
        if DaggerIter==1
            % set the expert's demonstration
            DCM.Ps = ExpertDemonstration;
            originPs = DCM.Ps;
        else

            for i=2:DaggerIter
                DCM.Ps(HFiter-(i-2),:,:) = outputCell{DaggerIter-(i-2),HFiter-(i-1)};
                disp(['HFiter:' num2str(HFiter) '  DaggerIter:' num2str(DaggerIter)]);
                disp(['DCM ' num2str(HFiter-(i-1)) ' outputcell:' num2str(DaggerIter-(i-2)) ' ' num2str(HFiter-(i-1)) ]);
            end
            
            originPs = DCM.Ps;
            expMask = rand(size(ExpertDemonstration,2),size(ExpertDemonstration,3))<expertRatio;
            
            for i=2:HFiter
                DCM.Ps(i,expMask) = ExpertDemonstration(i,expMask);
            end
            
        end

        disp(['put p from:' num2str(DaggerIter) ' ' num2str(HFiter-1)]);
        
          
        % Add new objective function into allFits
        allFits.addFitme(DCM);
       
        lb = 0*ones(size(Cs));
        ub = 1*ones(size(Cs));
        
        options = optimset('Jacobian','on','TolFun',1.0e-5, ...
           'TolX',3.0e-8,'MaxFunEvals',nIteration,'Display','iter-detailed'); 

        % Set the range of the objective function we want in allFits
        currentInd = size(allFits.fitmes,2);
        range = zeros(1,size(allFits.fitmes,2));
        range(startInd:currentInd) = 1;
        allFits.nSelect = range;
        
        disp(['DaggerIter:' num2str(DaggerIter) '  HFiter:' num2str(HFiter) '  start:' num2str(startInd)]);
       % disp(['currentInd:' num2str(currentInd) '  input:' currentInputStr '->' num2str(currentObjective)]);
        
        Cs
        %%%%%%%%%%%%%%%%%%% Optimization algorithm here
        [Cs,resnorm,residual,exitflag,output,lambda,jacobian] = ...
           lsqnonlin(@allFits.errJac,Cs,lb,ub,options);

        %allPars(:,HFiter) = pt;
        HFiter
        allCs{HFiter} = Cs;
        %lastX = pt;
        Cs

        
        DCM.Ps = originPs;
        
        [err,Jac] = DCM.errJac(Cs);
        % Store the output density matrices. These will be the input of next HFiter
        %outputCell{DaggerIter, HFiter} = allFits.fitmes{currentInd}.getAllRho();
        Ps_iter = DCM.currentPs;
        outputCell{DaggerIter+1, HFiter} = Ps_iter;
        
        disp(['store p to:' num2str(DaggerIter+1) ' ' num2str(HFiter)]);
      
    end

end


% Base line (0.5  0.5)

DCM2 = DistributeCoeff(Path,FragRangeTest,envs);
[err,Jac] = DCM2.errJac(1);
Cs = [0 1];
for j=2:HF_Iteration%/step_size
    % Now to get the list of latest Ps which were claculated
    Ps_iter = DCM2.currentPs;
    % This is a cell array of {nfrag,nenv+1}
    % Now add the Ps from iteration 1 to the Ps 
    DCM2.Ps(j,:,:) = Ps_iter;
    [err,Jac] = DCM2.errJac(Cs);
    Cs = zeros(1,j+1);
    Cs(end-1:end) = [0.5 0.5]
    BASE_LINE_ERROR(j) = sum(err.^2);
end

BASE_LINE_ERROR


% Evaluation (DAgger)
DCM3 = DistributeCoeff(Path,FragRangeTest,envs);
[err,Jac] = DCM3.errJac(1);

for j=2:HF_Iteration%/step_size
    % Now to get the list of latest Ps which were claculated
    Ps_iter = DCM3.currentPs;
    % This is a cell array of {nfrag,nenv+1}
    % Now add the Ps from iteration 1 to the Ps 
    DCM3.Ps(j,:,:) = Ps_iter;
    
    if(j>numel(allCs))
        Cs = zeros(1,j);
        Cs(end-1:end) = [0.5 0.5];
        [err,Jac] = DCM3.errJac(Cs);
    else
        [err,Jac] = DCM3.errJac(allCs{j});
    end
    
    
    ERROR(j) = sum(err.^2);
end

ERROR


% Original line (Pure Linear Combination)

DCM4 = DistributeCoeff(Path,FragRangeTest,envs);
[err,Jac] = DCM4.errJac(1);

for j=2:HF_Iteration%/step_size
    % Now to get the list of latest Ps which were claculated
    Ps_iter = DCM4.currentPs;
    % This is a cell array of {nfrag,nenv+1}
    % Now add the Ps from iteration 1 to the Ps 
    DCM4.Ps(j,:,:) = Ps_iter;
    
    
    if(j>numel(allCs))
        Cs = zeros(1,j);
        Cs(end-1:end) = [0.5 0.5];
        [err,Jac] = DCM4.errJac(Cs);
    else
        [err,Jac] = DCM4.errJac(oriCs{j});
    end
    
    
    Original_ERROR(j) = sum(err.^2);
end

disp('Naive 0.5 0.5')
BASE_LINE_ERROR

disp('DAgger with step size of 2')
ERROR

disp('Pure Linear Combination')
Original_ERROR

save('DAggerWithExpert_step2.mat','ERROR','Original_ERROR','BASE_LINE_ERROR','allCs', 'HF_Iteration','envs','Path','FragRange','FragRangeTest')


