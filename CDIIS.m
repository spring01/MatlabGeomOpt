classdef CDIIS < handle
    
    properties %(Access = private)
        
        fockVectors;
        errorCommutatorVectors;
        
        S_Half;
        
        startError = 0.1;
        
        densVectors;
        
    end
    
    methods
        
        function obj = CDIIS(overlapMatrix, numVectors)
            if(nargin < 2)
                numVectors = 5;
            end
            obj.fockVectors = zeros(numel(overlapMatrix), numVectors);
            obj.errorCommutatorVectors = zeros(numel(overlapMatrix), numVectors);
            
            obj.S_Half = sqrtm(overlapMatrix);
            
            % 
            obj.densVectors = zeros(numel(overlapMatrix), numVectors);
            %
        end
        
        function Push(obj, newFockVector, newDensVector)
            % push in matrices or rows or columns
            newFockVector = reshape(newFockVector, [], 1);
            newDensVector = reshape(newDensVector, [], 1);
            
            % push new Fock in
            obj.fockVectors(:, 1:end-1) = obj.fockVectors(:, 2:end);
            obj.fockVectors(:, end) = newFockVector;
            
            %
            obj.densVectors(:, 1:end-1) = obj.densVectors(:, 2:end);
            obj.densVectors(:, end) = newDensVector;
            %
            
            % push new commutator error in
            obj.errorCommutatorVectors(:, 1:end-1) = obj.errorCommutatorVectors(:, 2:end);
            
            FtDt = obj.S_Half ...
                \ reshape(newFockVector, sqrt(length(newFockVector)), []) ...
                * reshape(newDensVector, sqrt(length(newDensVector)), []) ...
                * obj.S_Half;
            obj.errorCommutatorVectors(:, end) = ...
                reshape(FtDt - FtDt', [], 1);
        end
        
        function better = IAmBetter(obj) % than EDIIS
            better = 0;
            if(sum(abs(obj.fockVectors(:,1))) ...
                    && obj.startError > max(abs(obj.errorCommutatorVectors(:, end))))
                better = 1;
            end
        end
        
        function newFockVector = Extrapolate(obj)
            onesVec = ones(obj.NumVectors(),1);
            hessian = [ ...
                obj.errorCommutatorVectors'*obj.errorCommutatorVectors, onesVec; ...
                onesVec', 0];
            if(rcond(hessian(1:end-1,1:end-1)) > 1e-15 && rcond(hessian) > 1e-15)
                diisCoefficients = hessian \ [zeros(obj.NumVectors(),1); 1];
                newFockVector = obj.fockVectors ...
                    * diisCoefficients(1:end-1);
                
%                 disp(diisCoefficients(1:end-1));
            else
                newFockVector = obj.fockVectors(:,end);
            end
        end
        
        %
        function newDensVector = ExtrapolateDensity(obj)
            onesVec = ones(obj.NumVectors(),1);
            hessian = [ ...
                obj.errorCommutatorVectors'*obj.errorCommutatorVectors, onesVec; ...
                onesVec', 0];
            if(rcond(hessian(1:end-1,1:end-1)) > 1e-15 && rcond(hessian) > 1e-15)
                diisCoefficients = hessian \ [zeros(obj.NumVectors(),1); 1];
                newDensVector = obj.densVectors ...
                    * diisCoefficients(1:end-1);
            else
                newDensVector = obj.densVectors(:,end);
            end
            
%             diisCoefficients = hessian \ [zeros(length(useIndices),1); 1];
%             newDensVector = obj.densVectors(:,useIndices) ...
%                 * diisCoefficients(1:end-1);
        end
        %
        
    end
    
    methods (Access = private)
        
        function num = NumVectors(obj)
            num = size(obj.fockVectors, 2);
        end
        
    end
    
end