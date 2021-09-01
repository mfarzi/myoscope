function thisP = singleObjMap(obj, thisX, p)
    % Helper function to map input x to p for a scalar LINKER object.
    switch obj.type
        case 'cos'
            lBound = obj.lowerBound;
            if obj.linked == 0
                uBound = obj.upperBound;
            elseif obj.linked == 1
                uBound = p(obj.crossLinkNo);
            elseif obj.linked == 2
                uBound = obj.upperBound - p(obj.crossLinkNo);
            else
                error('Unrecognised value %d for linked property', obj.linked);
            end
            thisP = cos(thisX)^2*(uBound-lBound)+lBound;

        case 'squared'
            thisP = thisX^2;

        case 'linear'
            thisP = thisX;

        case 'dummy'
            % note for dummy varialbes, thisX is an alias name for
            % thisP
            thisP = obj.fval(p, obj.crossLinkNo);

        case 'constant'
            thisP = obj.storedParamVal;

         otherwise
            error('The linking type %s is not recognised.', obj.type);
    end % switch obj.type
end % of scalarLink method