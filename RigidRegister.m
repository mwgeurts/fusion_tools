function rigid = RigidRegister(ref, new, varargin)
% RigidRegister rigidly registers a new CT image to a reference image. 
% This function supports both MATLAB and Plastimatch-based registration.
% The type of registration and basic parameters can be provided 
%
% If CT to density tables are provided, the reference image is converted 
% to new-equivalent by first interpolating to density using the reference 
% IVDT and then subsequently interpolating back to image value using the 
% new IVDT.
%
% The following variables are required for proper execution: 
%
%   ref:        structure containing the image data, dimensions, width,
%               start coordinates, and IVDT. See tomo_extract\LoadImage for 
%               an example of this format
%   new:        structure containing the image data, dimensions, width,
%               start coordinates and IVDT. See tomo_extract\LoadDailyImage
%               for an example of this format.
%
% In addition, the following options may be provided as name/value pairs.
% If not defined, the registration will be MATLAB using MSE.
%
%   method:     string contianing the algorithm to use for registration.  
%               Can be 'PLASTIMATCH' or 'MATLAB'
%   metric:     string containing metric. Can be 'MSE' 'GM' (plastimatch 
%               only) or 'MI'. If not provided, will default to 'MSE'
%   bone:       logical indicating whether to mask registration to only 
%               bony anatomy (true) or full image (false). If not provided,
%               will default to full image. This flag only has an effect if
%               IVDT tables are provided.
%   levels:     the number of levels to use during optimization.
%               Plastimatch is configured to support up to 3 levels.
%   iterations: number of iterations to run at each level.
%
% The following variables are returned upon succesful completion:
%
%   rigid:      six element vector of registration parameters [pitch yaw 
%               roll x y z], where angles are in degrees and distances are
%               in cm.
%
% Copyright (C) 2017 University of Wisconsin Board of Regents
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the  
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see http://www.gnu.org/licenses/.

% Initialize default options
method = 'MATLAB';
metric = 'MSE';
bone = false;
levels = 3;
iterations = 30;
bonethresh = 1.1;
bodythresh = 0.6;

% Loop through remaining input arguments
for i = 1:2:length(varargin)

    if strcmpi(varargin{i}, 'method')
        method = varargin{i+1};
    elseif strcmpi(varargin{i}, 'metric')
        metric = varargin{i+1};
    elseif strcmpi(varargin{i}, 'bone')
        bone = varargin{i+1};
    elseif strcmpi(varargin{i}, 'levels')
        levels = varargin{i+1};
    elseif strcmpi(varargin{i}, 'iterations')
        iterations = varargin{i+1};
    end
end

% If the bone flag is set
if exist('Event', 'file') == 2 && bone
    
    % Note use of bony anatomy to event log
    Event('Registering reference and new images using bony anatomy');
    
elseif exist('Event', 'file') == 2    
    
    % Note use of full image to event log
    Event('Registering reference and new images using full image');
end
    
% Log which registration method was chosen
if exist('Event', 'file') == 2  
    Event(['Method ', method, ' selected for registration']);
    
    % Start timer
    t = tic;
end

% If IVDT data exists
if isfield(ref, 'ivdt') && isfield(new, 'ivdt')
    
    % Convert reference image to equivalent new-IVDT image
    ref.data = interp1(new.ivdt(:,2), new.ivdt(:,1), interp1(ref.ivdt(:,1), ...
        ref.ivdt(:,2), ref.data, 'linear', 'extrap'), 'linear', 'extrap');

    % Note conversion in log
    if exist('Event', 'file') == 2  
        Event('Reference converted to equivalent image values using IVDT');
    end
end

% Execute registration based on method variable
switch method

%% Use Plastimatch 6-DOF rigid registration
% This rigid registration technique requires plastimatch be installed
% on this workstation.  Data is passed to/from plastimatch using the ITK 
% .mha file format and text command file. See 
% http://iopscience.iop.org/0031-9155/55/21/001 for additional details
% on the algorithm.
case 'PLASTIMATCH'
    
    %% Build reference MHA file
    % Generate a temprary filename for the reference image
    refFilename = [tempname, '.mha'];
    
    % Open a write file handle to the temporary reference image
    fid = fopen(refFilename, 'w', 'l');
    
    % Start writing the ITK header
    fprintf(fid, 'ObjectType=Image\n');
    fprintf(fid, 'NDims=3\n');
    
    % Specify the dimensions of the reference image
    fprintf(fid, 'DimSize=%i %i %i\n', ref.dimensions);
    
    % Specify the data format (USHORT referring to unsigned 16-bit integer)
    fprintf(fid,'ElementType=MET_USHORT\n');
    
    % Specify the byte order as little
    fprintf(fid,'ElementByteOrderMSB=False\n');
    
    % Specify the reference voxel widths (in mm)
    fprintf(fid, 'ElementSize=%i %i %i\n', ref.width*10);
    
    % Specify the reference voxel spacing to equal the widths (in mm)
    fprintf(fid, 'ElementSpacing=%i %i %i\n', ref.width*10);
    
    % Specify the coordinate frame origin (in mm)
    fprintf(fid, 'Origin=%i %i %i\n', ref.start*10);
    
    % Complete the .mha file header
    fprintf(fid, 'ElementDataFile=LOCAL\n');
    
    % Write the reference image data to the temporary file as uint16
    fwrite(fid, ref.data, 'ushort', 0, 'l');
    
    % Close the file handle
    fclose(fid);
    
    % Clear the temporary variable
    clear fid;
    
    % Log where the reference file was saved
    if exist('Event', 'file') == 2  
        Event(['Reference image written to ', refFilename]);
    end
    
    %% Build mask for reference image (excluding outside MVCT FOV)
    % Initialize null array of the same size as the reference image
    refMask = zeros(ref.dimensions);
    
    % Create meshgrid the same size as one image
    [x,y] = meshgrid(ref.start(1):ref.width(1):...
        ref.start(1) + ref.width(1) * ...
        (ref.dimensions(1) - 1), ref.start(2):...
        ref.width(2):ref.start(2) + ...
        ref.width(2) * (ref.dimensions(2) - 1));
    
    % Loop through each reference image slice
    for i = 1:ref.dimensions(3)
        
        % If the reference slice IEC-Y coordinate value is within the new
        % image slice range
        if ref.start(3)+(i*ref.width(3)) > ...
                new.start(3) && ref.start(3) + (i * ...
                ref.width(3)) < new.start(3) + ...
                new.dimensions(3) * new.width(3)
            
            % Set the mask to 1 within the new image FOV
            refMask(:,:,i) = sqrt(x.^2+y.^2) < newFOV/2 - 0.1;
        end
    end
    
    % If the bone flag is enabled, only include bone densities
    if isfield(ref, 'ivdt') && bone
        
        % Convert to density
        fixed = interp1(ref.ivdt(:,1), ref.ivdt(:,2), ref.data, ...
            'linear', 'extrap');
        
        % Update the fixed image to only include values above 1.10 g/cc
        refMask = refMask .* ceil((fixed - bonethresh) / 1e3);
        
    % Otherwise, include all densities above the threshold
    elseif isfield(ref, 'ivdt')
        
        % Convert to density
        fixed = interp1(ref.ivdt(:,1), ref.ivdt(:,2), ref.data, ...
            'linear', 'extrap');
        
        % Update the fixed image to only include values above 0.6 g/cc
        refMask = refMask .* ceil((fixed - bodythresh) / 1e3);
    end
    
    % Generate a temporary file name for the reference image mask
    refMaskFilename = [tempname, '.mha'];
    
    % Open a write file handle to the temporary reference image mask
    fid = fopen(refMaskFilename, 'w', 'l');
    
    % Start writing the ITK header
    fprintf(fid,'ObjectType=Image\n');
    fprintf(fid,'NDims=3\n');
    
    % Specify the dimensions of the reference image
    fprintf(fid, 'DimSize=%i %i %i\n', ref.dimensions);
    
    % Specify the data format (USHORT referring to unsigned 16-bit integer)
    fprintf(fid,'ElementType=MET_USHORT\n');
    
    % Specify the byte order as little
    fprintf(fid,'ElementByteOrderMSB=False\n');
    
    % Specify the reference voxel widths (in mm)
    fprintf(fid, 'ElementSize=%i %i %i\n', ref.width*10);
    
    % Specify the merged voxel spacing to equal the widths (in mm)
    fprintf(fid, 'ElementSpacing=%i %i %i\n', ref.width*10);
    
    % Specify the coordinate frame origin (in mm)
    fprintf(fid, 'Origin=%i %i %i\n', ref.start*10);
    
    % Complete the .mha file header
    fprintf(fid, 'ElementDataFile=LOCAL\n');
    
    % Write the reference image mask data to the temporary file as uint16
    fwrite(fid, refMask, 'ushort', 0, 'l');
    fclose(fid);
    Event(['Reference mask image written to ', refMaskFilename]); 
    
    %% Build new MHA file
    % Generate a temporary file name for the new image
    newFilename = [tempname, '.mha'];
    
    % Open a write file handle to the temporary new image
    fid = fopen(newFilename, 'w', 'l');
    
    % Start writing the ITK header
    fprintf(fid, 'ObjectType=Image\n');
    fprintf(fid, 'NDims=3\n');
    
    % Specify the dimensions of the merged image
    fprintf(fid, 'DimSize=%i %i %i\n', new.dimensions);
    
    % Specify the data format (USHORT referring to unsigned 16-bit integer)
    fprintf(fid, 'ElementType=MET_USHORT\n');
    
    % Specify the byte order as little
    fprintf(fid, 'ElementByteOrderMSB=False\n');
    
    % Specify the new voxel widths (in mm)
    fprintf(fid, 'ElementSize=%i %i %i\n', new.width*10);
    
    % Specify the new voxel spacing to equal the widths (in mm)
    fprintf(fid, 'ElementSpacing=%i %i %i\n', new.width*10);
    
    % Specify the coordinate frame origin (in mm)
    fprintf(fid, 'Origin=%i %i %i\n', new.start*10);
    
    % Complete the .mha file header
    fprintf(fid, 'ElementDataFile=LOCAL\n');
    
    % Write the merged image data to the temporary file as uint16
    fwrite(fid, new.data, 'uint16', 0, 'l');
    
    % Close the file handle
    fclose(fid);
    
    % Clear the temporary variable
    clear fid;
    
    % Log where the new file was saved
    Event(['New image written to ', newFilename]);
    
    %% Build mask for new image (excluding outside FOV)
    % Initialize null array of the same size as the new image
    newMask = zeros(new.dimensions);
    
    % Create meshgrid the same size as one image
    [x,y] = meshgrid(new.start(1):new.width(1):...
        new.start(1) + new.width(1) * ...
        (new.dimensions(1) - 1), new.start(2):...
        new.width(2):new.start(2) + ...
        new.width(2) * (new.dimensions(2) - 1));
    
    % Set the first mask slice to one within the FOV
    newMask(:,:,1) = sqrt(x.^2+y.^2) < newFOV/2 - 0.1;
    
    % Loop through each slice
    for i = 2:new.dimensions(3)
        
        % Copy the new mask to each slice
        newMask(:,:,i) = newMask(:,:,1);
    end
    
    % If the bone flag is enabled, only include bone densities
    if isfield(new, 'ivdt') && bone
        
        % Convert to density
        moving = interp1(new.ivdt(:,1), new.ivdt(:,2), new.data, ...
            'linear', 'extrap');
        
        % Update the fixed image to only include values above 1.10 g/cc
        newMask = newMask .* ceil((moving - bonethresh) / 1e3);
        
    % Otherwise, include all densities above the threshold
    elseif isfield(new, 'ivdt')
        
        % Convert to density
        moving = interp1(new.ivdt(:,1), new.ivdt(:,2), new.data, ...
            'linear', 'extrap');
        
        % Update the fixed image to only include values above 0.6 g/cc
        newMask = newMask .* ceil((moving - bodythresh) / 1e3);
    end
    
    % Generate a temporary file name for the new image mask
    newMaskFilename = [tempname, '.mha'];
    
    % Open a write file handle to the temporary new image mask
    fid = fopen(newMaskFilename, 'w', 'l');
    
    % Start writing the ITK header
    fprintf(fid, 'ObjectType=Image\n');
    fprintf(fid, 'NDims=3\n');
    
    % Specify the dimensions of the new image mask
    fprintf(fid, 'DimSize=%i %i %i\n', new.dimensions);
    
    % Specify the data format (USHORT referring to unsigned 16-bit integer)
    fprintf(fid, 'ElementType=MET_USHORT\n');
    
    % Specify the byte order as little
    fprintf(fid,'ElementByteOrderMSB=False\n');
    
    % Specify the new voxel widths (in mm)
    fprintf(fid, 'ElementSize=%i %i %i\n', new.width*10);
    
    % Specify the new voxel spacing to equal the widths (in mm)
    fprintf(fid, 'ElementSpacing=%i %i %i\n', new.width*10);
    
    % Specify the coordinate frame origin (in mm)
    fprintf(fid, 'Origin=%i %i %i\n', new.start*10);
    
    % Complete the .mha file header
    fprintf(fid, 'ElementDataFile=LOCAL\n');
    
    % Write the merged image mask data to the temporary file as uint16
    fwrite(fid, newMask, 'ushort', 0, 'l');
    
    % Close the file handle
    fclose(fid);
    
    % Clear the temporary variable
    clear fid;
    
    % Log where the new mask file was saved
    if exist('Event', 'file') == 2  
        Event(['New mask image written to ', newMaskFilename]);
    end
    
    %% Build plastimatch command file
    % Generate a temporary file name for the command file
    commandFile = [tempname, '.txt'];
    
    % Open a write file handle to the temporary command file
    fid = fopen(commandFile, 'w');
    
    % Specify the inputs to the registration
    fprintf(fid, '[GLOBAL]\n');
    fprintf(fid, 'fixed=%s\n', refFilename);
    fprintf(fid, 'moving=%s\n', newFilename);
    fprintf(fid, 'fixed_mask=%s\n', refMaskFilename);
    fprintf(fid, 'moving_mask=%s\n', newMaskFilename);
    
    % Generate a temporary filename for the resulting coefficients
    adjustments = [tempname, '.txt'];
    
    % Specify the output file
    fprintf(fid, 'xform_out=%s\n', adjustments);
    
    % If at least two levels are requested
    if levels >= 3
    
        % Specify stage 1 deformable image registration parameters.  Refer to 
        % http://plastimatch.org/registration_command_file_reference.html for
        % more information on these parameters
        fprintf(fid, '[STAGE]\n');
        fprintf(fid, 'xform=align_center\n');
    end
    
    % If at least two levels are requested
    if levels >= 2
        
        % Specify stage 2 parameters
        fprintf(fid, '[STAGE]\n');
        fprintf(fid, 'impl=plastimatch\n');
        fprintf(fid, 'xform=rigid\n');
        fprintf(fid, 'optim=versor\n');
        switch metric
            case 'MI'
                fprintf(fid, 'metric=mi\n');
            otherwise
                fprintf(fid, 'metric=mse\n');
        end
        fprintf(fid, 'max_its=%i\n', iterations);
        fprintf(fid, 'min_step=0.1\n');
        fprintf(fid, 'res=4 4 2\n');
        fprintf(fid, 'threading=cuda\n');
    end

    % Specify stage 3 parameters
    fprintf(fid, '[STAGE]\n');
    fprintf(fid, 'impl=plastimatch\n');
    fprintf(fid, 'xform=rigid\n');
    fprintf(fid, 'optim=versor\n');
    switch metric
        case 'MI'
            fprintf(fid, 'metric=mi\n');
        case 'GM'
            fprintf(fid, 'metric=gm\n');
        otherwise
            fprintf(fid, 'metric=mse\n');
    end
    fprintf(fid, 'max_its=%i\n', iterations);
    fprintf(fid, 'min_step=0.1\n');
    fprintf(fid, 'res=1 1 1\n');
    fprintf(fid, 'threading=cuda\n');
    
    %% Run plastimatch
    % Log execution of system call
    if exist('Event', 'file') == 2  
        Event(['Executing plastimatch register ', commandFile]);
    end
    
    % Execute plastimatch using system call, saving the output and status
    [status, cmdout] = system(['plastimatch register ', commandFile]);
    
    % If the status == 0, the command completed successfully
    if status == 0
        
        % Log output
        Event(cmdout);
    else
        % Otherwise, plastimatch didn't complete succesfully, so log the 
        % resulting command output as an error
        Event(cmdout, 'ERROR');
    end
    
    % Clear temporary variables
    clear status cmdout commandFile;
    
    %% Read in registration result
    % Open file handle to temporary file
    fid = fopen(adj, 'r');
    
    % Retrieve the first line of the result text
    tline = fgetl(fid);
    
    % Initialize temporary variables to flag if the results are found
    flag1 = 0;
    flag2 = 0;
    
    % Start a while loop to read in result text
    while ischar(tline)
        
        % Search for the text line containing the rigid registration,
        % storing the results and flag if the results are found
        [rigid, flag1] = sscanf(tline, ...
            '\nParameters: %f %f %f %f %f %f\n');
        
        % Search for the text line containing the rigid registration
        % origin, storing the origin and flag if the origin is found
        [origin, flag2] = sscanf(tline, '\nFixedParameters: %f %f %f\n');
        
        % Read in the next line of the results file
        tline = fgetl(fid);
    end
    
    % Close the file handle
    fclose(fid);
    
    % Clear the file handle
    clear fid;
    
    % If both flags are set, the results were successfully found
    if flag1 > 0 && flag2 > 0
        
        % Log an error indicating the the results were not parsed
        % correctly.  This usually indicates the registration failed
        if exist('Event', 'file') == 2  
            Event(['Unable to parse plastimatch results from ', adj], ...
                'ERROR'); 
        else
            error(['Unable to parse plastimatch results from ', adj]);
        end
        
    elseif exist('Event', 'file') == 2 
        
        % Otherwise, log success
        Event(['Plastimatch results read from ', adj]);
    end
    
    % If the registration origin is not equal to the DICOM center
    if ~isequal(origin, [0 0 0])
        
        % Log an error
        if exist('Event', 'file') == 2 
            Event(['Error: non-zero centers of rotation are not supported', ...
                ' at this time'], 'ERROR'); 
        else
            error(['Error: non-zero centers of rotation are not supported', ...
                ' at this time']);
        end
    end
    
    % Clear temporary variables
    clear flag1 flag2 origin;
    
    % Permute registration back to IEC coordinates
    rigid = -rigid([1 2 3 4 6 5]);
    
    % Convert rotations to degrees
    rigid(1:3) = rad2deg(rigid(1:3));
    
    % Report registration adjustments
    if exist('Event', 'file') == 2 
        Event(sprintf(['Rigid registration matrix [pitch yaw roll x y z] ', ...
            'computed as [%E %E %E %E %E %E] in %0.3f seconds'], ...
            rigid, toc(t)));
    end
    
    % Clear temporary variables
    clear referenceFilename newFilename referenceMaskFilename ...
        newMaskFilename adj commandFile;
    
%% Use MATLAB 6-DOF rigid registration
% This rigid registration technique requires MATLAB's Image Processing
% Toolbox imregtform function.
case 'MATLAB'
    
    % If the bone flag is enabled, only include bone densities
    if isfield(ref, 'ivdt') && bone
        
        % Convert to density
        fixed = interp1(ref.ivdt(:,1), ref.ivdt(:,2), ref.data, ...
            'linear', 'extrap');
        
        % Update the fixed image to only include values above 1.10 g/cc
        fixed = ref.data .* ceil((fixed - bonethresh) / 1e3);
        
    % Otherwise, include all densities above the threshold
    elseif isfield(ref, 'ivdt')
        
        % Convert to density
        fixed = interp1(ref.ivdt(:,1), ref.ivdt(:,2), ref.data, ...
            'linear', 'extrap');
        
        % Update the fixed image to only include values above 0.6 g/cc
        fixed = ref.data .* ceil((fixed - bodythresh) / 1e3);
    end
    
    %% Set fixed (reference) image and reference coordinates
    % Generate a reference meshgrid in the x, y, and z dimensions using the
    % start and width structure fields
    Rfixed = imref3d(size(ref.data), [ref.start(2) ...
        ref.start(2) + ref.width(2) * ...
        (size(ref.data, 2) - 1)], [ref.start(1) ...
        ref.start(1) + ref.width(1) * ...
        (size(ref.data, 1) - 1)], [ref.start(3) ...
        ref.start(3) + ref.width(3) * ...
        (size(ref.data, 3) - 1)]);
    
    %% Set moving (new) image and reference coordinates
    % If the bone flag is enabled, only include bone densities
    if isfield(new, 'ivdt') && bone
        
        % Convert to density
        moving = interp1(new.ivdt(:,1), new.ivdt(:,2), new.data, ...
            'linear', 'extrap');
        
        % Update the fixed image to only include values above 1.10 g/cc
        moving = new.data .* ceil((moving - bonethresh) / 1e3);
        
    % Otherwise, include all densities above the threshold
    elseif isfield(new, 'ivdt')
        
        % Convert to density
        moving = interp1(new.ivdt(:,1), new.ivdt(:,2), new.data, ...
            'linear', 'extrap');
        
        % Update the fixed image to only include values above 0.6 g/cc
        moving = new.data .* ceil((moving - bodythresh) / 1e3);
    end
    
    % Generate a reference meshgrid in the x, y, and z dimensions using the
    % start and width structure fields
    Rmoving = imref3d(size(new.data), [new.start(2) new.start(2) + ...
        new.width(2) * (size(new.data, 2) - 1)], [new.start(1) ...
        new.start(1) + new.width(1) * (size(new.data, 1) - 1)], ...
        [new.start(3) new.start(3) + new.width(3) * ...
        (size(new.data, 3) - 1)]);
    
    %% Run rigid registration
    % Initialize Regular Step Gradient Descent MATLAB object
    optimizer = registration.optimizer.RegularStepGradientDescent();
    
    % Set number of iterations to run
    optimizer.MaximumIterations = iterations;
    
    % Set maximum step size
    optimizer.MaximumStepLength = ref.width(1)/2;
    
    % Set minimum step size
    optimizer.MinimumStepLength = ref.width(1)/1e9;
    
    % Set metric
    switch metric
        
        % If the user chose Mutual Information
        case 'MI'
            
            % Initialize Mattes Mutual Information metric MATLAB object
            metric = registration.metric.MattesMutualInformation;

            % Log start of optimization
            if exist('Event', 'file') == 2 
                Event(['Executing imregtform rigid using Mattes ', ...
                    'mutual information']);
            end
            
        % Default to Mean Square Error
        otherwise
            
            % Initialize Mattes Mutual Information metric MATLAB object
            metric = registration.metric.MeanSquares;

            % Log start of optimization
            if exist('Event', 'file') == 2 
                Event('Executing imregtform rigid using mean square error');
            end
    end
    
    % Execute imregtform using 3 resampling levels
    tform = imregtform(moving, Rmoving, fixed, Rfixed, 'rigid', optimizer, ...
        metric, 'DisplayOptimization', 1, 'PyramidLevels', levels);
    
    % Clear temporary variables
    clear moving fixed Rmoving Rfixed metric optimizer;
    
    % Verify resulting transformation matrix is valid (the values (1,1) and
    % (3,3) must not be zero for atan2 to compute correctly)
    if tform.T(1,1) ~= 0 || tform.T(3,3) ~= 0

        % Compute pitch
        rigid(1) = -atan2(-tform.T(1,3), ...
            sqrt(tform.T(2,3)^2 + tform.T(3,3)^2));
                
        % Compute yaw
        rigid(2) = atan2(tform.T(2,3), tform.T(3,3));
        
        % Compute roll
        rigid(3) = -atan2(tform.T(1,2), tform.T(1,1));
    else
        % Otherwise, atan2 cannot compute, so throw an error
        if exist('Event', 'file') == 2 
            Event('Incompatible registration matrix determined', ...
                'ERROR');
        else
            error('Incompatible registration matrix determined');
        end
    end
    
    % Convert to degrees
    rigid(1:3) = rad2deg(rigid(1:3));
    
    % Set x, y, and z values
    rigid(4) = tform.T(4,2);
    rigid(5) = -tform.T(4,3);
    rigid(6) = -tform.T(4,1);
    
    % Clear transformation array
    clear tform;
    
    % Report registration results
    if exist('Event', 'file') == 2 
        Event(sprintf(['Rigid registration matrix [pitch yaw roll x y z] ', ...
            'computed as [%E %E %E %E %E %E] in %0.3f seconds'], ...
            rigid, toc(t)));
    end

% Otherwise, the method passed to RegisterImages was not supported    
otherwise
    
    % Throw error
    if exist('Event', 'file') == 2 
        Event(['Unsupported method ', method, ' passed to RegisterImages'], ...
            'ERROR');
    else
        error(['Unsupported method ', method, ' passed to RegisterImages']);
    end
end

% Clear temporary variables
clear t ref new method metric bone levels iterations bodythresh bonethresh;