classdef CFDIIS < handle
    
    properties (Access = private)
        
        fockVectors;
        errorVectors;
        
        S_Half;
        
        startError = 0.1;
        
        densVectors;
        finalDensVec;
        finalDensVecOrtho;
        
        weight;
        
    end
    
    methods
        
        function obj = CFDIIS(overlapMatrix, finalDensVector, weight, numVectors)
            if(nargin < 4)
                numVectors = 5;
            end
            obj.fockVectors = zeros(numel(overlapMatrix), numVectors);
            obj.errorVectors = zeros(2*numel(overlapMatrix), numVectors);
            
            obj.S_Half = sqrtm(overlapMatrix);
            
            % 
            obj.densVectors = zeros(numel(overlapMatrix), numVectors);
            obj.finalDensVec = reshape(finalDensVector, [], 1);
            obj.finalDensVecOrtho = reshape(obj.S_Half * reshape(finalDensVector, sqrt(numel(finalDensVector)), []) * obj.S_Half, [], 1);
            obj.weight = weight;
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
            obj.errorVectors(:, 1:end-1) = obj.errorVectors(:, 2:end);
            
            FtDt = obj.S_Half ...
                \ reshape(newFockVector, sqrt(length(newFockVector)), []) ...
                * reshape(newDensVector, sqrt(length(newDensVector)), []) ...
                * obj.S_Half;
            Dt = obj.S_Half * reshape(newDensVector, sqrt(length(newDensVector)), []) * obj.S_Half;
            
            
            obj.errorVectors(:, end) = ...
                [reshape(FtDt - FtDt', [], 1); obj.weight .* (reshape(Dt, [], 1) - obj.finalDensVecOrtho)];
%             obj.errorVectors(:, end) = ...
%                 [obj.weight .* (reshape(Dt, [], 1) - obj.finalDensVecOrtho)];
        end
        
        function better = IAmBetter(obj) % than EDIIS
            better = 0;
            if(sum(abs(obj.fockVectors(:,1))) ...
                    && obj.startError > max(abs(obj.errorVectors(:, end))))
                better = 1;
            end
        end
        
        %
        function newDensVector = ExtrapolateDensity(obj)
            onesVec = ones(obj.NumVectors(),1);
            hessian = [ ...
                obj.errorVectors'*obj.errorVectors, onesVec; ...
                onesVec', 0];
            useIndices = 1:obj.NumVectors();
            for i = 1:obj.NumVectors()-1
                if(rcond(hessian(1:end-1,1:end-1)) > 1e-15 && rcond(hessian) > 1e-15)
                    diisCoefficients = hessian \ [zeros(length(useIndices),1); 1];
                    newDensVector = obj.densVectors(:,useIndices) ...
                        * diisCoefficients(1:end-1);
                    return;
                else
                    hessian = hessian(2:end, 2:end);
                    useIndices = useIndices(2:end);
                end
            end
            newDensVector = obj.densVectors(:,end);
        end
        %
        
    end
    
    methods (Access = private)
        
        function num = NumVectors(obj)
            num = size(obj.fockVectors, 2);
        end
        
    end
    
end