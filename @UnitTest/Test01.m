function Test01(testCase)
% UNIT TEST 1: MATLAB rigid registration accuracy
%
% DESCRIPTION: This unit test runs a series of random registration
% adjustments against RigidRegister using the MATLAB algorithm. Using a
% randomly generated registration vector, each reference image is
% transformed using RigidTransform. The original and transformed images are 
% then passed to RigidRegister, and the resulting registration vector 
% compared to the known modification. In essence, this test verifies that 
% RigidRegister and RigidTransform are self-consistent.
%
% RELEVANT REQUIREMENTS: None
%
% INPUT DATA: DICOM test image dataset
%
% CONDITION A (+): The RMS error should be less than 0.1 cm and 0.1 deg for
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
        
        % Execute RigidTransform
        mod = RigidTransform(image, v);
        
        % Execute RigidRegister
        r = RigidRegister(image, mod, 'method', 'MATLAB');
        
        % Store errors
        errors = vertcat(errors, reshape((r + v)', 3, 2)); %#ok<AGROW>
    end
end

% Verify the results are not empty
testCase.verifyNotEmpty(errors);

% Compute RMS error
rmsdeg = sqrt(mean(errors(:,1).^2));
rmscm = sqrt(mean(errors(:,2).^2));

% Verify results are less than tolerance
testCase.verifyLessThan(rmsdeg, 0.1);
testCase.verifyLessThan(rmscm, 0.1);

% Store results
testCase.StoreResults('results', sprintf('%0.2f deg<br>%0.2f cm',...
    rmsdeg, rmscm));