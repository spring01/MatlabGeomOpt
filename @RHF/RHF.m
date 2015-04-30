classdef RHF < handle
    
    properties (SetAccess = protected)
        
        overlapMat;
        coreHamilt;
        nucRepEnergy;
        numElectrons;
        matpsi2;
        
        hfEnergy;
        densMat;
        
    end
    
    properties (Access = protected)
        
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
        
        [hfEnergy, iter] = SCF(obj, iniDensMat)
        
    end
    
    methods (Access = protected)
        
        function gMat = DensToG(obj, densMat)
            gMat = 2 .* obj.matpsi2.JK_DensToJ(densMat) ... % +2J
                - obj.matpsi2.JK_DensToK(densMat); % -K
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