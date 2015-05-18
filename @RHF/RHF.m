classdef RHF < handle
    
    properties (SetAccess = protected)
        
        orbital;
        densVec;
        
        finalFockVec;
        
    end
    
    properties (Access = protected)
        
        overlapMat;
        coreHamilt;
        nucRepEnergy;
        numElectrons;
        matpsi2;
        
        maxSCFIter = 500;
        RMSDensityThreshold = 1e-8;
        MaxDensityThreshold = 1e-6;
        EnergyThreshold = 1e-6;
        
    end
    
    methods
        
        function obj = RHF(properties)
            obj.overlapMat = properties.overlapMat;
            obj.coreHamilt = properties.coreHamilt;
            obj.matpsi2 = properties.matpsi2;
            obj.nucRepEnergy = properties.nucRepEnergy;
            obj.numElectrons = properties.numElectrons;
        end
        
    end
    
    methods (Access = protected)
        
        function gMat = DensToG(obj, densMat)
            gMat = 2 .* obj.matpsi2.JK_DensToJ(densMat) ... % +2J
                - obj.matpsi2.JK_DensToK(densMat); % -K
        end
        
        function [densVec, elecEnergy, orbital, orbEigValues] ...
                = DiagonalizeFock(obj, fockMat, inv_S_Half)
            [orbitalOtho, orbEigValues] = eig(inv_S_Half*fockMat*inv_S_Half);
            [orbEigValues, ascend_order] = sort(diag(orbEigValues));
            orbital = inv_S_Half * orbitalOtho(:, ascend_order);
            occOrb = orbital(:, 1:obj.numElectrons/2);
            densVec = reshape(occOrb * occOrb', [], 1);
            elecEnergy = sum(orbEigValues(1:obj.numElectrons/2));
        end
        
    end
    
    methods (Static)
                
        function rhf = MatPsi2Interface(matpsi2)
            properties.overlapMat = matpsi2.Integrals_Overlap();
            properties.coreHamilt = matpsi2.Integrals_Kinetic() + matpsi2.Integrals_Potential();
            properties.matpsi2 = matpsi2;
            properties.nucRepEnergy = matpsi2.Molecule_NucRepEnergy();
            properties.numElectrons = matpsi2.Molecule_NumElectrons();
            rhf = RHF(properties);
        end
        
    end
    
end