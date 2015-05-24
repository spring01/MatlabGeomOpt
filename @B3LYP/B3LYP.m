classdef B3LYP < handle
    
    properties (SetAccess = protected)
        
        orbital;
        densVec;
        
        finalFockVec;
        
        energySet;
        
    end
    
    properties (Access = protected)
        
        overlapMat;
        coreHamilt;
        nucRepEnergy;
        numElectrons;
        matpsi2;
        
        maxSCFIter = 500;
        RMSDensityThreshold = 1e-8;
        MaxDensityThreshold = 1e-1;
        EnergyThreshold = 1e-6;
        
        currentV;
        
    end
    
    methods
        
        function obj = B3LYP(properties)
            obj.overlapMat = properties.overlapMat;
            obj.coreHamilt = properties.coreHamilt;
            obj.matpsi2 = properties.matpsi2;
            obj.matpsi2.DFT_Initialize('b3lyp');
            obj.nucRepEnergy = properties.nucRepEnergy;
            obj.numElectrons = properties.numElectrons;
        end
        
    end
    
    methods (Access = protected)
        
        function gMat = DensToG(obj, dens)
            obj.currentV = obj.matpsi2.DFT_DensToV(dens);
            gMat = 2 .* obj.matpsi2.JK_DensToJ(dens) ...
                - 0.2 .* obj.matpsi2.JK_DensToK(dens) ...
                + obj.currentV;
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
                
        function obj = MatPsi2Interface(matpsi2)
            properties.overlapMat = matpsi2.Integrals_Overlap();
            properties.coreHamilt = matpsi2.Integrals_Kinetic() + matpsi2.Integrals_Potential();
            properties.matpsi2 = matpsi2;
            properties.nucRepEnergy = matpsi2.Molecule_NucRepEnergy();
            properties.numElectrons = matpsi2.Molecule_NumElectrons();
            obj = B3LYP(properties);
        end
        
    end
    
end