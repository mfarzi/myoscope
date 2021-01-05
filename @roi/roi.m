classdef roi < handle
    % ROI
    %
    %   A ROI object facilitate the management and visualisation of voxels
    %   in a 3D or 2D image. 
    %
    %   properties: (private)
    %       nrow              - number of rows in the image
    %       ncol              - number of columns in the image
    %       ndep              - number of stacks in the third dimension
    %       index             - index to the selected voxels in the ROI 
    %
    %   methods (public):
    %       getIndex          - get index 
    %       getVoxel          - return coordinates for voxels in the ROI
    %       update            - update ROI voxels
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (Access = 'private')                              
        index = [];           % index to the selected voxels in ROI
                              % i + nrow*(j-1) + nrow*ncol*(k-1)     
        nrow = [];            % number of rows in the image
        ncol = [];            % number of columns in the image
        ndep = [];            % number of voxels along the third direction 
    end
    
    methods 
        function obj = roi(varargin)
            %ROI Construct Function.
            %
            %   roi(bw) construct an ROI object with input binary mask
            %
            %   roi(nr,nc,nd) construct a 3D ROI object of size [nr,nc,nd]
            %   where all voxels are selected as true
            %
            %   roi(nr,nc) construct a 2D ROI object of size [nr,nc]
            %   where all voxels are selected as true
            %           
            % Parse the parameter-value pairs
            if nargin==0
                % do nothing. Default Construction Function
            else
                [idx, nrow, ncol, ndep] = obj.parseInput(varargin{:});
                obj.nrow = nrow;
                obj.ncol = ncol;
                obj.ndep = ndep;
                obj.index = idx;
            end
        end % of constructor
        
        function ind = getIndex(obj)
            % return index
            ind = obj.index;
        end
        
        function sub = getSubscript(obj)
            if obj.ndep==1 % 2D mask
                [x, y] = ind2sub([obj.nrow, obj.ncol], obj.index);
                sub = [x;y];
            else %3D mask
                [x, y, z] = ind2sub([obj.nrow, obj.ncol, obj.ndep], obj.index);
                sub = [x; y; z];
            end
        end
        
        function d = getSize(obj)
            if obj.ndep==1 %2D mask
                d = [obj.nrow, obj.ncol];
            else %3D mask
                d = [obj.nrow, obj.ncol, obj.ndep];
            end
        end
        
        function write(obj, filename)
            % save roi as a text file
            [filepath,name,ext] = fileparts(filename);
            if isempty(filepath)
                filepath = pwd;
            end
            if ~strcmp(ext, '.roi')
                warning('MATLAB:ROI:save',...
                        'The filename extension is changed to ".roi".');
                ext = '.roi';    
            end
            filename = fullfile(filepath, strcat(name,ext));
            if isfile(filename)
                warning('MATLAB:ROI:save',...
                        'Filename is overwritten.\n %s', filename);
            end
            fileId = fopen(filename,'w');
            fprintf(fileId, '##ROI\n');
            fprintf(fileId, '#$nrow\n');
            fprintf(fileId, '%d\n', obj.nrow);
            fprintf(fileId, '#$ncol\n');
            fprintf(fileId, '%d\n', obj.ncol);
            fprintf(fileId, '#$ndep\n');
            fprintf(fileId, '%d\n', obj.ndep);
            fprintf(fileId, '#$index\n');
            fprintf(fileId, '%d ', obj.index);
            fclose(fileId);
        end%of method save
                
        function [C, iA, iB] = intersect(obj, thisRoi)
            if ~isa(thisRoi, 'roi')
                error('MTLAB:ROI:intersect',...
                      'Input should be of type "roi"');
            end
            assert(thisRoi.nrow==obj.nrow && ...
                   thisRoi.ncol==obj.ncol && ...
                   thisRoi.ndep==obj.ndep, 'MATLAB:ROI:update', ...
                   'The size of input roi does not match.');
            [idx, iA, iB] = intersect(obj.index, thisRoi.index);
            C = roi(obj.nrow, obj.ncol, obj.ndep);
            C.setIndex(idx);
        end
        
        function [C, iA, iB] = union(obj, thisRoi)
            if ~isa(thisRoi, 'roi')
                error('MTLAB:ROI:intersect',...
                      'Input should be of type "roi"');
            end
            assert(thisRoi.nrow==obj.nrow && ...
                   thisRoi.ncol==obj.ncol && ...
                   thisRoi.ndep==obj.ndep, 'MATLAB:ROI:update', ...
                   'The size of input roi does not match.');
            [idx, iA, iB] = union(obj.index, thisRoi.index);
            C = roi(obj.nrow, obj.ncol, obj.ndep);
            C.setIndex(idx);
        end
        
        function bw = getMask(obj)
            if obj.ndep==1 %2D mask
                bw = false(obj.nrow, obj.ncol);
            else
                bw = false(obj.nrow, obj.ncol, obj.ndep);
            end
            bw(obj.index) = true;
        end
        
        function setIndex(obj, idx)
            % check idx is row or column vector
            if isrow(idx)
                obj.index = idx;
            elseif iscolumn(idx)
                obj.index = idx';
            else
                error('MATLAB:ROI:setIndex',...
                      'Indput idex shuld be a row or column vector.');
            end
        end%of setIndex
            
        function setSubscript(obj, sub)
            % check subscript is a 2D matrix with two or three rows
            assert(ismatrix(sub), 'MATLAB:ROI:setSubscript',...
                    'Input subscript should be a 2D matrix.');
            if obj.ndep==1 %2D image
                assert(size(sub,1)==2, 'MATLAB:ROI:setSubscript',...
                    'Input subscripts should have two rows.');
            else
                assert(size(sub,1)==3, 'MATLAB:ROI:setSubscript',...
                    'Input subscripts should have three rows.');
            end

            if obj.ndep==1 %2D image
                obj.index = sub2ind([obj.nrow, obj.ncol],sub(1,:),sub(2,:));
            else
                obj.index = sub2ind([obj.nrow, obj.ncol, obj.ndep],...
                                sub(1,:), sub(2,:), sub(3,:));
            end
        end%of setSubscript
        %\\
    end % of methods (public)
    
    methods (Static)
        function obj = read(filename)
            obj = roi();
            [filepath,name,ext] = fileparts(filename);
            if isempty(filepath)
                filepath = pwd;
                filename = fullfile(filepath, strcat(name, ext));
            end
            if ~isfile(filename)
                error('MATLAB:ROI:load',...
                      'Filename does not exist.\n%s', filename);
            end
            fileId = fopen(filename, 'r');
            thisLine = fgetl(fileId);
            while ischar(thisLine)
                if strcmp(thisLine, '#$nrow')
                    obj.nrow = str2double(fgetl(fileId));
                end
                
                if strcmp(thisLine, '#$ncol')
                    obj.ncol = sscanf(fgetl(fileId),'%d');
                end
                
                if strcmp(thisLine, '#$ndep')
                    obj.ndep = sscanf(fgetl(fileId),'%d');
                end
                
                if strcmp(thisLine, '#$index')
                    obj.index = sscanf(fgetl(fileId),'%d')';
                end
                
                thisLine = fgetl(fileId);
            end
            fclose(fileId);
        end%of method read
    end
        
    
    methods (Static, Access = 'private')
        function [idx, nrow, ncol, ndep] = parseInput(varargin)
            switch nargin
                case 1
                    if isa(varargin{1}, 'char')
                        [idx, nrow, ncol, ndep] = read(varargin{1});
                    else
                        bw = varargin{1};
                        assert(islogical(bw), 'MATLAB:ROI:parseInput', ...
                               strcat('Class ROI with %s',...
                               ' input is not identified.'), class(bw));
                        nrow = size(bw, 1);
                        ncol = size(bw, 2);
                        ndep = size(bw, 3);
                        idx = find(bw)';
                    end
                case 2
                    nrow = varargin{1};
                    assert(isscalar(nrow) && isnumeric(nrow) && rem(nrow, 1)==0,...
                           'MATLAB:ROI:parseInput', ...
                           strcat('roi(nrow, ncol, ndep) is defined',...
                           ' for scalar integer inputs only.'));
                       
                    ncol = varargin{2};
                    assert(isscalar(ncol) && isnumeric(ncol) && rem(ncol, 1)==0,...
                           'MATLAB:ROI:parseInput', ...
                           strcat('roi(nrow, ncol, ndep) is defined',...
                           ' for scalar integer inputs only.'));
                       
                    ndep = 1;
                    
                    idx = 1:nrow*ncol;
                case 3
                    nrow = varargin{1};
                    assert(isscalar(nrow) && isnumeric(nrow) && rem(nrow, 1)==0,...
                           'MATLAB:ROI:parseInput', ...
                           strcat('roi(nrow, ncol, ndep) is defined',...
                           ' for scalar integer inputs only.'));
                       
                    ncol = varargin{2};
                    assert(isscalar(ncol) && isnumeric(ncol) && rem(ncol, 1)==0,...
                           'MATLAB:ROI:parseInput', ...
                           strcat('roi(nrow, ncol, ndep) is defined',...
                           ' for scalar integer inputs only.'));
                       
                    ndep = varargin{3};
                    assert(isscalar(ndep) && isnumeric(ndep) && rem(ndep, 1)==0,...
                           'MATLAB:ROI:parseInput', ...
                           strcat('roi(nrow, ncol, ndep) is defined', ...
                           ' for scalar integer inputs only.'));

                    idx = 1:nrow*ncol*ndep;
                otherwise
                    error('MATLAB:ROI:parseInput',...
                          strcat('Class ROI with this input argument',...
                          ' is not identified. Use roi(bw) or',...
                          ' roi(nrow, ncol) or roi(nrow, ncol, ndep).'));
            end%of switch case
        end%of parseInput    
        %\\
    end % of methods (static, private)
    %\\
end%of class roi