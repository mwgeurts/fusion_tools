function image = LoadDICOMImages(path)
% LoadDICOMImages scans a provided directory for enclosed DICOM files. 
% DICOM files are loaded and returned as structures. See the dicom_tools
% submodule functions for more information on the structure format.
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

% Scan the directory for DICOM files
if exist('Event', 'file') == 2
    Event(['Scanning ', path, ' for DICOM files'], 'UNIT');
end

% Retrieve folder contents of selected directory
list = dir(path);

% Initialize folder counter
i = 0;

% Initialize empty 3D array for images and vector of slice locations
% (the data may not be loaded in correct order; these will be used to
% re-sort the slices later)
images = [];
sliceLocations = [];

% Start recursive loop through each folder, subfolder
while i < length(list)

    % Increment current folder being analyzed
    i = i + 1;

    % If the folder content is . or .., skip to next folder in list
    if strcmp(list(i).name, '.') || strcmp(list(i).name, '..')
        continue

    % Otherwise, if the folder content is a subfolder    
    elseif list(i).isdir == 1

        % Retrieve the subfolder contents
        sublist = dir(fullfile(path, list(i).name));

        % Look through the subfolder contents
        for j = 1:size(sublist, 1)

            % If the subfolder content is . or .., skip to next subfolder 
            if strcmp(sublist(j).name, '.') || ...
                    strcmp(sublist(j).name, '..')
                continue
            else

                % Otherwise, replace the subfolder name with its full
                % reference
                sublist(j).name = fullfile(list(i).name, ...
                    sublist(j).name);
            end
        end

        % Append the subfolder contents to the main folder list
        list = vertcat(list, sublist); %#ok<AGROW>

        % Clear temporary variable
        clear sublist;

    % Otherwise, see if the file is a DICOM file
    else

        % Attempt to parse the DICOM header
        try
            % Execute dicominfo
            info = dicominfo(fullfile(path, list(i).name));

            % Verify storage class field exists
            if ~isfield(info, 'MediaStorageSOPClassUID')
                continue
            end

            % If CT or MR, add to imagefiles
            if strcmp(info.MediaStorageSOPClassUID, ...
                    '1.2.840.10008.5.1.4.1.1.2') || ...
                    strcmp(info.MediaStorageSOPClassUID, ...
                    '1.2.840.10008.5.1.4.1.1.4')
                
                % Append this slice's location to the sliceLocations vector
                sliceLocations(length(sliceLocations)+1) = ...
                    info.ImagePositionPatient(3); %#ok<*AGROW>

                % Append this slice's image data to the images array
                images(size(images,1)+1,:,:) = dicomread(info); %#ok<*AGROW>
            end

        % If an exception occurs, the file is not a DICOM file so skip
        catch
            continue
        end
    end
end

% Initialize empty return argument
image = struct();

% Only continue if at least one image was found
if isempty(images)
    return
end

% Retrieve start voxel IEC-X coordinate from DICOM header, in cm
image.start(1) = info.ImagePositionPatient(1) / 10 * ...
    info.ImageOrientationPatient(1);

% Adjust IEC-Z to inverted value, in cm
image.start(2) = -(info.ImagePositionPatient(2) * ...
    info.ImageOrientationPatient(5) + info.PixelSpacing(2) * ...
    (size(images, 2) - 1)) / 10; 

% Retrieve X/Z voxel widths from DICOM header, in cm
image.width(1) = info.PixelSpacing(1) / 10;
image.width(2) = info.PixelSpacing(2) / 10;

% If patient is Head First
if isequal(info.ImageOrientationPatient, [1;0;0;0;1;0]) || ...
        isequal(info.ImageOrientationPatient, [-1;0;0;0;-1;0]) 

    % Sort sliceLocations vector in descending order
    [~, indices] = sort(sliceLocations, 'descend');

    % Store start voxel IEC-Y coordinate, in cm
    image.start(3) = -max(sliceLocations) / 10;

% Otherwise, if the patient is Feet First
elseif isequal(info.ImageOrientationPatient, [-1;0;0;0;1;0]) || ...
        isequal(info.ImageOrientationPatient, [1;0;0;0;-1;0]) 

    % % Sort sliceLocations vector in ascending order
    [~,indices] = sort(sliceLocations, 'ascend');

    % Store start voxel IEC-Y coordinate, in cm
    image.start(3) = min(sliceLocations) / 10;
end 

% Compute slice location differences
widths = diff(sliceLocations(indices));

% Store mean slice position difference as IEC-Y width, in cm
image.width(3) = abs(mean(widths)) / 10;

% Initialize daily image data array as single type
image.data = single(zeros(size(images, 3), size(images, 2), ...
    size(images, 1)));

% Loop through each slice
for i = 1:length(sliceLocations)

    % Set the image data based on the index value
    image.data(:, :, i) = ...
        single(rot90(permute(images(indices(i), :, :), [2 3 1])));
end

% Remove values less than zero (some DICOM images place voxels outside the
% field of view to negative values)
image.data = max(image.data, 0);

% Flip images in IEC-X direction
image.data = flip(image.data, 1);

% Create dimensions structure field based on the daily image size
image.dimensions = size(image.data);