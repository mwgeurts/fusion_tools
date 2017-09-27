function mod = RigidTransform(image, v)
% RigidTransform transforms the provided image by a rigid registration
% vector, returning a modified image that has the same coordinate system
% and voxel positions as the original image. Image values outside of the
% original image are set to zero.
%
% The following variables are required for proper execution: 
%
%   image:  structure containing the image data, dimensions, width,
%           start coordinates, and IVDT.  See tomo_extract or dicom_tools
%           submodules for more documentation on this file format.
%   v:      6 element registration vector in [pitch yaw roll x y z] IEC
%           coordinates, or a 4x4 tranformation matrix.
%
% The following variables are returned upon succesful completion:
%
%   mod:    structure containing a transformed image
%           (converted back to the daily IVDT), transformation matrix, 
%           dimensions, width, and start coordinates
%
% Author: Mark Geurts, mark.w.geurts@gmail.com
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

%% Determine registration
% If the variable input is a registration vector
if isvector(v) && length(v) == 6
    
    %% Generate transformation matrix
    % Log start of transformation
    if exist('Event', 'file') == 2        
        Event('Generating transformation matrix');
    end
    
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
    
% Otherwise, if input is a transformation matrix
elseif ismatrix(v) && isequal(size(v), [4 4])
    
    % Store input
    tform = v;
    
% Otherwise, throw an error
else
    if exist('Event', 'file') == 2
        Event('Incorrect second argument. See documentation for formats', ...
            'ERROR');
    else
        error('Incorrect second argument. See documentation for formats');
    end
end

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
mod = image;

% Store transformation matrix
mod.tform = tform;

% Use try-catch statement to attempt to perform interpolation using GPU.  
% If a GPU compatible device is not available (or fails due to memory), 
% automatically revert to CPU based technique
try
    % Initialize and clear GPU memory
    gpuDevice(1);

    % Interpolate the daily image dataset to the reference dataset's
    % transformed coordinates using GPU linear interpolation, and store to 
    % varargout{1}
    mod.data = gather(interp3(gpuArray(secX), gpuArray(secY), ...
        gpuArray(secZ), gpuArray(image.data), gpuArray(refX), ...
        gpuArray(refY), gpuArray(refZ), 'linear', 0));

    % Clear memory
    gpuDevice(1);

catch

    % Interpolate the daily image dataset to the reference dataset's
    % transformed coordinates using CPU linear interpolation, and store to 
    % varargout{1}
    mod.data = interp3(refX, refY, refZ, image.data, secX, ...
        secY, secZ, '*linear', 0);
end

% Clear temporary variables
clear result ref1 refX refY refZ secX secY secZ;
