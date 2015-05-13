classdef CheatDIIS < handle
    
    properties %(Access = private)
        
        finalDensVector;
        densVectors;
        
    end
    
    methods
        
        function obj = CheatDIIS(finalDensity, numVectors)
            if(nargin < 2)
                numVectors = 5;
            end
            obj.finalDensVector = reshape(finalDensity, [], 1);
            obj.densVectors = zeros(numel(finalDensity), numVectors);
        end
        
        function Push(obj, newDensVector)
            % push in matrices or rows or columns
            newDensVector = reshape(newDensVector, [], 1);
            
            % push newDensVector in
            obj.densVectors(:,1:end-1) = obj.densVectors(:,2:end);
            obj.densVectors(:,end) = newDensVector;
        end
        
        function [diisCoefficients, predDensVector] = Predict(obj)
            onesVec = ones(obj.NumVectors(), 1);
            hessian = [obj.densVectors'*obj.densVectors, onesVec; ...
                onesVec', 0];
            useIndices = 1:obj.NumVectors();
            for i = 1:obj.NumVectors()-1
                if(rcond(hessian) > 1e-15)
                    break;
                else
                    hessian = hessian(2:end, 2:end);
                    useIndices = useIndices(2:end);
                end
            end
            diisCoefficients = hessian \ [obj.densVectors(:,useIndices)'*obj.finalDensVector; 1];
            predDensVector = obj.densVectors(:,useIndices) ...
                * diisCoefficients(1:end-1);
            diisCoefficients = diisCoefficients(1:end-1);
        end
        
%         function [diisCoefficients, predDensVector] = Predict(obj)
%             onesVec = ones(obj.NumVectors(), 1);
%             hessian = [obj.densVectors'*obj.densVectors];
%             useIndices = 1:obj.NumVectors();
%             for i = 1:obj.NumVectors()-1
%                 if(rcond(hessian) > 1e-15)
%                     break;
%                 else
%                     hessian = hessian(2:end, 2:end);
%                     useIndices = useIndices(2:end);
%                 end
%             end
%             diisCoefficients = hessian \ [obj.densVectors(:,useIndices)'*obj.finalDensVector];
%             predDensVector = obj.densVectors(:,useIndices) ...
%                 * diisCoefficients;
%         end
        
    end
    
    methods (Access = private)
        
        function num = NumVectors(obj)
            num = size(obj.densVectors, 2);
        end
        
    end
    
end
