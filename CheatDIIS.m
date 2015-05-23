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
        
%         function predDensVector = ExtrapolateDensity(obj)
%             onesVec = ones(obj.NumVectors(), 1);
%             hessian = [obj.densVectors'*obj.densVectors, onesVec; ...
%                 onesVec', 0];
%             useIndices = 1:obj.NumVectors();
%             for i = 1:obj.NumVectors()-1
%                 if(rcond(hessian) > 1e-15)
%                     break;
%                 else
%                     hessian = hessian(2:end, 2:end);
%                     useIndices = useIndices(2:end);
%                 end
%             end
%             diisCoefficients = hessian \ [obj.densVectors(:,useIndices)'*obj.finalDensVector; 1];
%             predDensVector = obj.densVectors(:,useIndices) ...
%                 * diisCoefficients(1:end-1);
%             diisCoefficients = diisCoefficients(1:end-1);
%             
%             disp(norm(predDensVector - obj.finalDensVector));
%         end
%         
%         function predDensVector = ExtrapolateDensity(obj)
%             hessian = obj.densVectors'*obj.densVectors;
%             useIndices = 1:obj.NumVectors();
%             for i = 1:obj.NumVectors()-1
%                 if(rcond(hessian) > 1e-15)
%                     diisCoefficients = hessian \ obj.densVectors(:,useIndices)'*obj.finalDensVector;
%                     predDensVector = obj.densVectors(:,useIndices) * diisCoefficients;
%                     disp(norm(predDensVector - obj.finalDensVector));
%                     return;
%                 else
%                     hessian = hessian(2:end, 2:end);
%                     useIndices = useIndices(2:end);
%                 end
%             end
%             predDensVector = obj.densVectors(:,end);
% 
%         end
        
        function predDensVector = ExtrapolateDensity(obj)
            XTrain = obj.densVectors;
            hessian = XTrain'*XTrain;
            useIndices = 1:obj.NumVectors();
            for iter = 1:size(XTrain, 2)-1
                if(rcond(hessian) > 1e-15)
                    break;
                else
                    useIndices = useIndices(2:end);
                end
            end
            
            XTrain = XTrain(:, useIndices);
            
            newXTrain = zeros(size(XTrain, 1), 0);
            for col1 = 1:size(XTrain, 2)
%                 for col2 = col1:min(col1+1, size(XTrain, 2))
                    newXTrain = [newXTrain, XTrain(:, col1).^2, XTrain(:, col1).^3, XTrain(:, col1).^4];
%                 end
            end
            XTrain = [XTrain, newXTrain];
            
            hessian = XTrain'*XTrain;
            diisCoefficients = hessian \ XTrain'*obj.finalDensVector;
            predDensVector = XTrain * diisCoefficients;
            disp(norm(predDensVector - obj.finalDensVector));
        end
        
    end
    
    methods (Access = private)
        
        function num = NumVectors(obj)
            num = size(obj.densVectors, 2);
        end
        
    end
    
end
