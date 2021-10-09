function [model, params, rmse, flag, sig, schemeFile] = runProcedure(filename, varargin)
    % runProcedure is a static method for the abstract class COMPARTMENT
    % 
    % runProcedure(filename) return the compartmen model, fitted
    % parameters, data, schemeFile, rmse, and exit flag during
    % optimisation.
    %       Input args:
    %         filename: Filename of type char to read the procedure.
    %
    %      Output Args:
    %            model: A compartment object
    %           params: A numeric matrix [nParams x nVoxels]
    %              sig: A numeric matrix [nScheme x nVoxels]
    %       schemeFile: A scheme object with nScheme measurements
    %             rmse: Root mean squared error
    %             flag: Exit flag for the optimisation
    %       
    %
    % see also: readProcedure, runProcedure
    
    [model, sig, schemeFile] = compartment.readProcedure(filename);
    
    nD = sig(1,1);
    roi = [nD; sig(2:nD+1,1); sig(nD+2,:)'];
    
    sig = sig(nD+3:end, :);
    
    [params, rmse, flag] = model.fit(sig, schemeFile, varargin{:});
    
    model.writeParams(filename, params, rmse, roi);
end