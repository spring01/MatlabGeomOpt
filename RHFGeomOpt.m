classdef RHFGeomOpt < handle
    
    properties (SetAccess = protected)
        
        initialMol;
        finalMol;
        
        rhf;
        
    end
    
    properties (Access = protected)
        
        matpsi2;
        thresMaxForce = 0.000450
        thresRMSForce = 0.000300;
        thresMaxDisplace = 0.001800;
        thresRMSDisplace = 0.001200;
        maxIter = 500;
        stepSize = 1;
        
    end
    
    methods
        
        function obj = RHFGeomOpt(molecule, basisName)
            obj.initialMol = molecule;
            obj.matpsi2 = MatPsi2(molecule.cartesian, basisName);
            obj.rhf = RHF.MatPsi2Interface(obj.matpsi2);
        end
        
        function [currGeom, iter] = DoGeomOpt(obj)
            [~, iterRHF] = obj.rhf.SCF();
            disp(iterRHF)
            currGeom = obj.matpsi2.Molecule_Geometry();
            obj.matpsi2.SCF_RunSCF();
            currForceVec = reshape(obj.matpsi2.SCF_Gradient(), [], 1);
            currHess = eye(length(currForceVec));
            for iter = 1:obj.maxIter
                maxForce = max(abs(currForceVec));
                rmsForce = sqrt(currForceVec' * currForceVec ./ length(currForceVec));
                
                nextGeom = currGeom - obj.stepSize .* reshape(currHess\currForceVec, [], 3);
                
                obj.matpsi2.Molecule_SetGeometry(nextGeom);
                
                prevOrbital = obj.rhf.orbital;
                obj.rhf = RHF.MatPsi2Interface(obj.matpsi2);
                [~, iterRHF] = obj.rhf.SCF(prevOrbital);
                obj.matpsi2.SCF_SetGuessOrb(prevOrbital);
                obj.matpsi2.SCF_RunSCF();
                nextForceVec = reshape(obj.matpsi2.SCF_Gradient(), [], 1);
                
                
                deltaForceVec = nextForceVec - currForceVec;
                deltaGeomVec = reshape(nextGeom - currGeom, [], 1);
                temp = deltaGeomVec'*currHess;
                nextHess = currHess ...
                    + (deltaForceVec*deltaForceVec')./(deltaForceVec'*deltaGeomVec) ...
                    - (temp'*temp)./(deltaGeomVec'*currHess*deltaGeomVec);
                
                displaceVec = reshape(nextGeom - currGeom, [], 1);
                maxDisplace = max(abs(displaceVec));
                rmsDisplace = sqrt(displaceVec' * displaceVec ./ length(displaceVec));
                if(maxForce < obj.thresMaxForce ...
                        && rmsForce < obj.thresRMSForce ...
                        && maxDisplace < obj.thresMaxDisplace ...
                        && rmsDisplace < obj.thresRMSDisplace)
                    obj.finalMol = Molecule([obj.initialMol.cartesian(:,1), nextGeom.*obj.initialMol.Bohr2Angstrom]);
                    break;
                end
                currGeom = nextGeom;
                currForceVec = nextForceVec;
                currHess = nextHess;
                
                disp(iterRHF)
                disp(maxForce)
                disp(rmsForce)
                disp(maxDisplace)
                disp(rmsDisplace)
            end
        end
        
    end
    
end
