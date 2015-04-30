classdef RHFGeomOpt < handle
    
    properties (SetAccess = protected)
        
        initialMol;
        finalMol;
        
    end
    
    properties %(Access = protected)
        
        matpsi2;
        thresMaxForce = 0.000450
        thresRMSForce = 0.000300;
        thresMaxDisplace = 0.001800;
        thresRMSDisplace = 0.001200;
        maxIter = 500;
        stepSize = 1.5;
        
    end
    
    methods
        
        function obj = RHFGeomOpt(molecule, basisName)
            obj.initialMol = molecule;
            obj.matpsi2 = MatPsi2(molecule.cartesian, basisName);
        end
        
        function [currGeom, iter] = DoGeomOpt(obj)
            currGeom = obj.matpsi2.Molecule_Geometry();
            obj.matpsi2.SCF_RunRHF();
            forceVec = reshape(obj.matpsi2.SCF_Gradient(), [], 1);
            currHess = diag(abs(forceVec));
            for iter = 1:obj.maxIter
                
                
                maxForce = max(abs(forceVec));
                rmsForce = sqrt(forceVec' * forceVec ./ length(forceVec));
                nextGeom = currGeom - obj.stepSize .* reshape(forceVec, [], 3);
                displaceVec = reshape(nextGeom - currGeom, [], 1);
                maxDisplace = max(abs(displaceVec));
                rmsDisplace = sqrt(displaceVec' * displaceVec ./ length(displaceVec));
                if(maxForce < obj.thresMaxForce ...
                        && rmsForce < obj.thresRMSForce ...
                        && maxDisplace < obj.thresMaxDisplace ...
                        && rmsDisplace < obj.thresRMSDisplace)
                    obj.finalMol = Molecule([obj.initialMol(:,1), currGeom.*obj.initialMol.Bohr2Angstrom]);
                    break;
                end
                obj.matpsi2.Molecule_SetGeometry(nextGeom);
                obj.matpsi2.SCF_RunRHF();
                forceVec = reshape(obj.matpsi2.SCF_Gradient(), [], 1);
                currGeom = nextGeom;
                disp(iter)
                disp(maxForce)
                disp(rmsForce)
                disp(maxDisplace)
                disp(rmsDisplace)
            end
        end
        
    end
    
end
