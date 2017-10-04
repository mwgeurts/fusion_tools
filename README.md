## Rigid and Deformable Registration Tools for MATLAB&reg;

by Mark Geurts <mark.w.geurts@gmail.com>
<br>Copyright &copy; 2017, University of Wisconsin Board of Regents

This package contains MATLAB&reg; wrapper functions for several [plastimatch<sup>&copy;</sup>](http://plastimatch.org) executables to register and manipulate volumetric images stored in the data format used by other modules, including [dicom_tools](https://github.com/mwgeurts/dicom_tools) and [tomo_extract](https://github.com/mwgeurts/tomo_extract). This is a submodule for other tools such as [mvct_retro](https://github.com/mwgeurts/mvct_retro). MATLAB is a registered trademark of MathWorks Inc. Plastimatch is Copyrighted by The General Hospital Corporation and is available at http://plastimatch.org.

At this time only rigid registration is available. Check back later as this module is actively being developed to add more functionality! 

## Contents

* [Installation and Use](README.md#installation-and-use)
* [Compatibility and Requirements](README.md#compatibility-and-requirements)
* [Tools and Examples](README.md#tools-and-examples)
  * [RigidRegister](README.md#rigidregister)
  * [RigidTransform](README.md#rigidtransform)
  * [FunctionWhich](README.md#functionwhich)
* [Event Calling](README.md#event-calling)
* [Unit Testing](README.md#unit-testing)
* [License](README.md#license)

## Installation and Use

To install the Rigid and Deformable Registration Tools, copy all MATLAB .m files and contents of the `bin64` folder from this repository into your MATLAB path. If installing as a submodule into another git repository, execute `git submodule add https://github.com/mwgeurts/fusion_tools`. 

The `RigidRegister()` function uses `FunctionWhich()` to determine which plastimatch executable to run, depending on system environment. Compiled versions of plastimatch 1.6.4 are provided for Windows 7, MacOS 10.10, and Ubuntu 14.04. To compile a new version of plastimatch, follow the [getting started instructions](http://plastimatch.org/getting_started.html) provided by the code owner and then replace the corresponding entry (Windows, MacOS, or Linux) in the `paths` variable in `FunctionWhich()`.

## Compatibility and Requirements

The Rigid and Deformable Registration Tools have been validated for MATLAB versions 9.1 through 9.2 on Macintosh OSX 10.10 (Yosemite),  Windows 7, and Ubuntu 14.04. These tools optionally use the Image Processing Toolbox MATLAB function `imregtform()` when performing MATLAB-based rigid registration and the Parallel Processing Toolbox `interp3()` when performing GPU-based interpolation. If this toolbox or a compatible GPU device is not found, the `RigidTransform()` function will automatically revert to CPU interpolation. Plastimatch registration does not require any additional toolboxes.

## Tools and Examples

The following subsections describe what inputs and return variables are used, and provides examples for basic operation of each tool. For more information, refer to the documentation within the source code.

### RigidRegister

`RigidRegister()` rigidly registers a target CT image to a reference image. This function supports both MATLAB and Plastimatch-based registration. The type of registration and basic parameters can be provided 

If CT to density tables are provided, the reference image is converted to target-equivalent by first interpolating to density using the reference IVDT and then subsequently interpolating back to image value using the target IVDT.

The following variables are required for proper execution: 

* ref: structure containing the image data, dimensions, width, start coordinates, and IVDT. See tomo_extract or dicom_tools submodules for more documentation on this file format.
* target: structure containing the image data, dimensions, width, start coordinates and IVDT. See tomo_extract or dicom_tools submodules for more documentation on this file format.

In addition, the following options may be provided as name/value pairs. If not defined, the registration will be MATLAB using MSE.

* method: string contianing the algorithm to use for registration. Can be 'PLASTIMATCH' or 'MATLAB'
* metric: string containing metric. Can be 'MSE', 'GM' (plastimatch only) or 'MI'. If not provided, will default to 'MSE'
* bone: logical indicating whether to mask registration to only bony anatomy (true) or full image (false). If not provided, will default to full image. This flag only has an effect if IVDT tables are provided.
* levels: the number of levels to use during optimization. For Plastimatch, each level multiplies the voxel size by 2
* iterations: number of iterations to run at each level.
* rotations: flag indicating whether to include rotations (true) or only optimize on translations (false)

The following variables are returned upon succesful completion:

* rigid: six element vector of registration parameters [pitch yaw roll x y z], where angles are in degrees and distances are in cm.

Below is an example of how this function is used. The [dicom_tools](https://github.com/mwgeurts/dicom_tools) submodule is used to load the images:

```matlab
% Add the dicom_tools functions
addpath('/path/to/dicom_tools/');

% Load DICOM images 
path = '/path/to/files/';  
names = { 
    '2.16.840.1.114362.1.5.1.0.101218.5981035325.299641582.274.1.dcm' 
    '2.16.840.1.114362.1.5.1.0.101218.5981035325.299641582.274.2.dcm'
    '2.16.840.1.114362.1.5.1.0.101218.5981035325.299641582.274.3.dcm'  
};  
imageA = LoadDICOMImages(path, names);  

path = '/path/to/files/';  
names = { 
    '2.16.840.1.114362.1.5.1.0.101218.5981035325.299641582.123.1.dcm' 
    '2.16.840.1.114362.1.5.1.0.101218.5981035325.299641582.123.2.dcm'
    '2.16.840.1.114362.1.5.1.0.101218.5981035325.299641582.123.3.dcm'  
};  
imageB = LoadDICOMImages(path, names); 

% Register the images using MATLAB MI-based rigid registration
rigid = RigidRegister(imageA, imageB, 'method', 'MATLAB', 'metric', 'MI');

% Re-register the images using PLASTIMATCH MSE-based rigid registration
rigid = RigidRegister(imageA, imageB, 'method', 'PLASTIMATCH');

% Re-register the images using only translations
rigid = RigidRegister(imageA, imageB, 'rotations', false);
```

### RigidTransform

`RigidTransform()` transforms the provided image by a rigid registration vector, returning a modified image that has the same coordinate system and voxel positions as the original image. Image values outside of the original image are set to zero.

This function optionally uses the Parallel Processing Toolbox `interp3()` when performing GPU-based interpolation. If this toolbox or a compatible GPU device is not found, it will automatically revert to CPU interpolation. 

The following variables are required for proper execution: 

* image: structure containing the image data, dimensions, width, start coordinates, and IVDT.  See tomo_extract or dicom_tools submodules for more documentation on this file format.
* v: 6 element registration vector in [pitch yaw roll x y z] IEC coordinates, or a 4x4 tranformation matrix.

The following variables are returned upon succesful completion:

* mod: structure containing a transformed image (converted back to the daily IVDT), transformation matrix, dimensions, width, and start coordinates

Below is an example of how this function is used alongside `RigidRegister()`. The [dicom_tools](https://github.com/mwgeurts/dicom_tools) submodule is used to load the images:

```matlab
% Add the dicom_tools functions
addpath('/path/to/dicom_tools/');

% Load DICOM images 
path = '/path/to/files/';  
names = { 
    '2.16.840.1.114362.1.5.1.0.101218.5981035325.299641582.274.1.dcm' 
    '2.16.840.1.114362.1.5.1.0.101218.5981035325.299641582.274.2.dcm'
    '2.16.840.1.114362.1.5.1.0.101218.5981035325.299641582.274.3.dcm'  
};  
imageA = LoadDICOMImages(path, names);  

path = '/path/to/files/';  
names = { 
    '2.16.840.1.114362.1.5.1.0.101218.5981035325.299641582.123.1.dcm' 
    '2.16.840.1.114362.1.5.1.0.101218.5981035325.299641582.123.2.dcm'
    '2.16.840.1.114362.1.5.1.0.101218.5981035325.299641582.123.3.dcm'  
};  
imageB = LoadDICOMImages(path, names); 

% Register the images using PLASTIMATCH MSE-based rigid registration
rigid = RigidRegister(imageA, imageB, 'method', 'PLASTIMATCH');

% Transform image B using the rigid registration results
modImageB = RigidTransform(imageB, rigid);
```

### FunctionWhich

`FunctionWhich()` attempts to find the location of the provided command. This function is part of the fusion_tools submodule, and is called by the other functions to determine how to execute plastimatch. The command is first searched for in the system path; if not found, the function will successively search through each path provided in the paths variable until it is found. This function is compatible with linux, MacOS, and Windows 7 and later operating systems. 

The following variables are required for proper execution: 

* command: string containing the command to be executed (plastimatch)

The following variables are returned upon successful completion:

* executable: string containing the path and command name to the command. If the function was not found, this string will be empty.

Below is an example of how this function is used:
 
```matlab
% Search for the function 'plastimatch'
executable = FunctionWhich('plastimatch');

% If found, use it to execute a registration
if ~isempty(executable) {
    system([executable, ' register commandfile']);
}
```

## Event Calling

These functions optionally return execution status and error information to an `Event()` function. If available in the MATLAB path, `Event()` will be called with one or two variables: the first variable is a string containing the status information, while the second is the status classification (WARN or ERROR). If the status information is only informative, the second argument is not included.  Finally, if no `Event()` function is available errors will still be thrown via the standard `error()` MATLAB function.

## Unit Testing

This function includes a `UnitTest` class to validate key features of the provided functions. When executed via `runtests()`, this class will automatically run a series of unit tests against the specified executable, and create Markdown reports summarizing the results. The MATLAB profiler is turned on during unit test execution to quantify the code coverage. Refer to the [MathWorks documentation](https://www.mathworks.com/help/matlab/class-based-unit-tests.html) for more information on running class-based unit tests.

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
