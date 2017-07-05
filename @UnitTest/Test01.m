function Test01(testCase)
% UNIT TEST 1: MATLAB rigid registration accuracy
%
% DESCRIPTION: This unit test runs a series of random registration
% adjustments against RigidRegister using the MATLAB algorithm. Using a
% randomly generated registration vector, each reference image is
% transformed using linear interpolation. The original and transformed
% images are then passed to RigidRegister, and the resulting registration
% vector compared to the known modification.
%
% RELEVANT REQUIREMENTS: None
%
% INPUT DATA: DICOM test datasets
%
% CONDITION A (+): The RMS error should be less than 0.5 cm and 0.5 deg for
% offsets and rotations, respectively.

% Log test
if exist('Event', 'file') == 2
    Event('Executing unit test 1', 'UNIT');
end

% Define number of random tests to run
n = 5;

% Initialize errors array
errors = [];

% Store test summary
testCase.StoreResults('summary', ...
    sprintf('MATLAB rigid registration RMS error (n = %i)', ...
    n * size(testCase.inputData, 1)));

% Loop through test datasets
for i = 1:size(testCase.inputData, 1)

    % Load DICOM image
    image = testCase.LoadDICOMImages(fullfile(testCase.inputData{i,2}));
    
    % Loop through each error
    for j = 1:n
    
        % Compute normally random registration matrix
        v = randn(1,6);
        
        % Generate 4x4 transformation matrix given a 6 element vector of 
        % [pitch yaw roll x y z].  
        tform(1,1) = cosd(v(3)) * cosd(v(1));
        tform(2,1) = cosd(v(3)) * sind(v(1)) * sind(v(2)) - sind(v(3)) * cosd(v(2));
        tform(3,1) = cosd(v(3)) * sind(v(1)) * cosd(v(2)) + sind(v(3)) * sind(v(2));
        tform(4,1) = v(6);
        tform(1,2) = sind(v(3)) * cosd(v(1));
        tform(2,2) = sind(v(3)) * sind(v(1)) * sind(v(2)) + cosd(v(3)) * cosd(v(2));
        tform(3,2) = sind(v(3)) * sind(v(1)) * cosd(v(2)) - cosd(v(3)) * sind(v(2));
        tform(4,2) = v(4);
        tform(1,3) = -sind(v(1));
        tform(2,3) = cosd(v(1)) * sind(v(2));
        tform(3,3) = cosd(v(1)) * cosd(v(2));
        tform(4,3) = -v(5);
        tform(1,4) = 0;
        tform(2,4) = 0;
        tform(3,4) = 0;
        tform(4,4) = 1;
        
        % Generate x, y, and z grids using start and width structure fields 
        % (which are stored in [x,z,y] format)
        [refX, refY, refZ] = meshgrid(image.start(2) + ...
            image.width(2) * (size(image.data, 2) - 1): ...
            -image.width(2):image.start(2), ...
            image.start(1):image.width(1):image.start(1) ...
            + image.width(1) * (size(image.data, 1) - 1), ...
            image.start(3):image.width(3):image.start(3) ...
            + image.width(3) * (size(image.data, 3) - 1));
        
        % Generate unity matrix of same size as reference data to aid in matrix 
        % transform
        ref1 = ones(image.dimensions);

        % Separately transform each reference x, y, z point by shaping all to 
        % vector form and dividing by transformation matrix
        result = [reshape(refX,[],1) reshape(refY,[],1) reshape(refZ,[],1) ...
            reshape(ref1,[],1)] / tform;

        % Reshape transformed x, y, and z coordinates back to 3D arrays
        secX = reshape(result(:,1), image.dimensions);
        secY = reshape(result(:,2), image.dimensions);
        secZ = reshape(result(:,3), image.dimensions);
        
        % Copy image
        modimage = image;
        
        % Use try-catch statement to attempt to perform interpolation using GPU.  
        % If a GPU compatible device is not available (or fails due to memory), 
        % automatically revert to CPU based technique
        try
            % Initialize and clear GPU memory
            gpuDevice(1);

            % Interpolate the daily image dataset to the reference dataset's
            % transformed coordinates using GPU linear interpolation, and store to 
            % varargout{1}
            modimage.data = gather(interp3(gpuArray(secX), gpuArray(secY), ...
                gpuArray(secZ), gpuArray(image.data), gpuArray(refX), ...
                gpuArray(refY), gpuArray(refZ), 'linear', 0));

            % Clear memory
            gpuDevice(1);

        catch
            
            % Interpolate the daily image dataset to the reference dataset's
            % transformed coordinates using CPU linear interpolation, and store to 
            % varargout{1}
            modimage.data = interp3(refX, refY, refZ, image.data, secX, ...
                secY, secZ, '*linear', 0);

        end

        % Clear temporary variables
        clear result ref1 refX refY refZ secX secY secZ;
        
        % Execute RigidRegister
        r = RigidRegister(image, modimage, 'method', 'MATLAB');
        
        % Store errors
        errors = vertcat(errors, reshape((r + v)', 3, 2)); %#ok<AGROW>
    end
end

% Store results
testCase.StoreResults('results', sprintf('%0.2f deg<br>%0.2f cm',...
    sqrt(mean(errors(:,1).^2)), sqrt(mean(errors(:,2).^2))));