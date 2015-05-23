classdef ComDIIS < handle
    
    properties %(Access = private)
        
        fockVectors;
        densVectors;
        S_Half;
        
    end
    
    methods
        
        function obj = ComDIIS(overlapMatrix, numVectors)
            if(nargin < 2)
                numVectors = 5;
            end
            obj.fockVectors = zeros(numel(overlapMatrix), numVectors);
            obj.densVectors = zeros(numel(overlapMatrix), numVectors);
            obj.S_Half = sqrtm(overlapMatrix);
        end
        
        function Push(obj, newFockVector, newDensVector)
            % push in matrices or rows or columns
            newFockVector = reshape(newFockVector, [], 1);
            newDensVector = reshape(newDensVector, [], 1);
            
            % push new Fock in
            obj.fockVectors(:, 1:end-1) = obj.fockVectors(:, 2:end);
            obj.fockVectors(:, end) = newFockVector;
            
            obj.densVectors(:, 1:end-1) = obj.densVectors(:, 2:end);
            obj.densVectors(:, end) = newDensVector;
        end
        
        function newDensVector = ExtrapolateDensity(obj)
            fockVecs = obj.fockVectors;
            densVecs = obj.densVectors;
            for i = 1:size(obj.fockVectors, 2)
                if(norm(densVecs(:, 1)) == 0)
                    newDensVector = obj.densVectors(:, end);
                    return;
                end
            end
            
            
            nbf = sqrt(size(fockVecs, 1));
            nVecs = size(fockVecs, 2);
            
            e = cell(nVecs, nVecs);
            for i = 1:nVecs
                for j = 1:nVecs
                    fock_i = reshape(fockVecs(:, i), nbf, []);
                    dens_j = reshape(densVecs(:, j), nbf, []);
                    FtDt = obj.S_Half \ (fock_i * dens_j) * obj.S_Half;
                    e{i, j} = FtDt - FtDt';
                end
            end
            
            H = zeros(nVecs, nVecs, nVecs, nVecs);
            
            
            for i = 1:nVecs
                for j = 1:nVecs
                    for k = 1:nVecs
                        for l = 1:nVecs
                            H(i,j,k,l) = sum(sum(e{i, j} .* e{k, l}));
                        end
                    end
                end
            end
            
            iniCoeffs = zeros(nVecs, 1);
            iniCoeffs(end) = 1;
            
            options = optimoptions(@fmincon, 'Display', 'off', 'GradObj', 'on', 'TolFun', 1e-14, 'TolCon', 1e-14);
            finalCoeffs = fmincon(@(coeffs)Target(coeffs, H), ...
                iniCoeffs, [], [], ones(1, length(iniCoeffs)), 1, [], [], [], options);
            disp(finalCoeffs)
            newDensVector = densVecs * finalCoeffs;
        end
        
    end
    
    methods (Access = private)
        
        function num = NumVectors(obj)
            num = size(obj.fockVectors, 2);
        end
        
    end
    
end