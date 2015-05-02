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
        stepSize = 1.5;
        
    end
    
    methods
        
        function obj = RHFGeomOpt(molecule, basisName)
            obj.initialMol = molecule;
            obj.matpsi2 = MatPsi2(molecule.cartesian, basisName);
            obj.rhf = RHF.MatPsi2Interface(obj.matpsi2);
        end
        
        function [currGeom, iter] = RunGeomOpt(obj)
            [~, iterRHF] = obj.rhf.SCF();
            disp(iterRHF)
            currGeom = obj.matpsi2.Molecule_Geometry();
            obj.matpsi2.SCF_SetSCFType('RHF');
            obj.matpsi2.SCF_RunSCF();
            currForceVec = reshape(obj.matpsi2.SCF_Gradient(), [], 1);
            currHess = eye(length(currForceVec));
            for iter = 1:obj.maxIter
                % next geometry
                displaceVec = obj.stepSize .* (-currHess\currForceVec);
                nextGeom = currGeom + reshape(displaceVec, [], 3);
                
                % SCF and gradient of next geometry
                obj.matpsi2.Molecule_SetGeometry(nextGeom);
                prevOrbital = obj.rhf.orbital;
                obj.rhf = RHF.MatPsi2Interface(obj.matpsi2);
                [~, iterRHF] = obj.rhf.SCF(prevOrbital);
                obj.matpsi2.SCF_SetGuessOrb(obj.rhf.orbital);
                obj.matpsi2.SCF_RunSCF();
                nextForceVec = reshape(obj.matpsi2.SCF_Gradient(), [], 1);
                
                % update Hessian
                deltaForceVec = nextForceVec - currForceVec;
                deltaGeomVec = reshape(nextGeom - currGeom, [], 1);
                temp = currHess * deltaGeomVec;
                currHess = currHess ...
                    + (deltaForceVec*deltaForceVec')./(deltaForceVec'*deltaGeomVec) ...
                    - (temp*temp')./(deltaGeomVec'*currHess*deltaGeomVec);
                
                % convergence criteria
                maxForce = max(abs(nextForceVec));
                rmsForce = sqrt(mean(nextForceVec.^2));
                maxDisplace = max(abs(displaceVec));
                rmsDisplace = sqrt(mean(displaceVec.^2));
                
                disp(iter)
                disp(iterRHF)
                disp(maxForce)
                disp(rmsForce)
                disp(maxDisplace)
                disp(rmsDisplace)
                
                if(maxForce < obj.thresMaxForce ...
                        && rmsForce < obj.thresRMSForce ...
                        && maxDisplace < obj.thresMaxDisplace ...
                        && rmsDisplace < obj.thresRMSDisplace)
                    obj.finalMol = Molecule([obj.initialMol.cartesian(:,1), nextGeom.*obj.initialMol.Bohr2Angstrom]);
                    break;
                end
                
                % start from next geometry
                currGeom = nextGeom;
                currForceVec = nextForceVec;
            end
        end
        
    end
    
end
