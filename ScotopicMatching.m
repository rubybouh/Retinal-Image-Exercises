% First need to create a rod spectral sensitivity function
% Wavelength vs. probability of absorption per unit of time 

% Get scotopic matching system matrix by putting in 1s at each entry
% and zeros in all the others
% Loop through all vectors

% And then this matrix should be a scaled version of your rod spectral
% sensitivity function

clear all
close all

%%% Creating rod spectral sensitivity function %%%
rng(0); % Set the random seed for reproducibility
wavelength = linspace(400, 700, 31);  % Wavelengths from 400 to 700 nm
rodSensitivity = rand(1, length(wavelength)); % Random power values for rod spectral sensitivity function

figure(1);
hold on
plot(wavelength, rodSensitivity);
xlabel('Wavelength (nm)');
ylabel('Relative Power');

%%% Finding the scotopic matching system matrix %%%

scotopicMatchingSystemMatrix = zeros(1, length(wavelength));
scalingFactor = 2; % Adjust as needed

for i = 1:length(wavelength) % Looping through wavelengths
    unitTestLight = zeros(length(wavelength), 1); 
    unitTestLight(i) = 1; % Set intensity of test light at the current wavelength to 1
    unitTestLight = unitTestLight .* scalingFactor; % Scale if desired
    systemMatrixEntry = rodSensitivity * unitTestLight;
    scotopicMatchingSystemMatrix(1, i) = systemMatrixEntry;
end

%%% Plotting scotopic matching system matrix %%%
% Should be the same as rod spectral sensitivity up to a linear
% transformation
plot(wavelength, scotopicMatchingSystemMatrix, 'LineWidth', 2, 'LineStyle', '--');
legend('Rod Spectral Sensitivity', 'Scotopic Matching System Matrix');
    

