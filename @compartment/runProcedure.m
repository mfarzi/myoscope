function [model, params, sig, schemeFile, rmse, flag] = runProcedure(filename)
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
    
    [params, rmse, flag] = model.fit(sig, schemeFile);
    model.writeParams(filename, params, rmse, flag);
end