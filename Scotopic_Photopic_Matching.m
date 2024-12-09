clear all
close all

%%%%% SCOTOPIC Wavelength Encoding %%%%%

%%% Creating rod spectral sensitivity function and primary light spectral distribution %%%
rng(0); % Set the random seed for reproducibility
wavelength = linspace(400, 700, 31);  % Wavelengths from 400 to 700 nm
rodSensitivity = rand(1, length(wavelength)); % Random power values for rod spectral sensitivity function
% Same as the system matrix of the rhodopsin absoprtion experiment
rng(15); % Set a different random seed 
primaryDist = rand(1, length(wavelength)) * 1/4; % Random power values for the primary light spectral distribution
figure(1);
plot(wavelength, primaryDist);
xlabel('Wavelength (nm)');
ylabel('Relative Power');
title('Primary (Match) Light Spectral Distribution')
figure(2);
hold on
plot(wavelength, rodSensitivity);
xlabel('Wavelength (nm)');
ylabel('Relative Power');

%%% Finding the scotopic matching system matrix %%%

scotopicMatchingSystemMatrix = zeros(1, length(wavelength));
scalingFactor = 2; % Adjust as needed
unscaledMatchIntensity = rodSensitivity * (primaryDist)';

for i = 1:length(wavelength) % Looping through wavelengths
    unitTestLight = zeros(length(wavelength), 1); 
    unitTestLight(i) = 1; % Set intensity of test light at the current wavelength to 1
    intensity = rodSensitivity * unitTestLight;
    knobSetting = intensity / unscaledMatchIntensity;
    scotopicMatchingSystemMatrix(1, i) = knobSetting;
end

%%% Plotting scotopic matching system matrix %%%
% Should be the same as rod spectral sensitivity up to a linear
% transformation
plot(wavelength, scotopicMatchingSystemMatrix, 'LineWidth', 2, 'LineStyle', '--');
legend('Rod Spectral Sensitivity', 'Scotopic Matching System Matrix');
    
%%%%% PHOTOPIC Wavelength Encoding %%%%%

%%% Creating cone spectral sensitivity function and primary light spectral distributions %%%
rng(0); % Set the random seed for reproducibility
wavelength = linspace(400, 700, 31);  % Wavelengths from 400 to 700 nm

L_cone = rand(1, length(wavelength)); % Random power values for cone spectral sensitivity function
M_cone = rand(1, length(wavelength));
S_cone = rand(1, length(wavelength));

coneSensitivity = [L_cone', M_cone', S_cone'];
% Synonym to "the system matrix of the cone photopigment absoprtion experiment"

rng(15); % Set a different random seed 
primaryDist1 = rand(1, length(wavelength)) * 1/4; % Random power values for the primary light spectral distribution
primaryDist2 = rand(1, length(wavelength)) * 1/4;
primaryDist3 = rand(1, length(wavelength)) * 1/4;

primaryDist = [primaryDist1', primaryDist2', primaryDist3'];

figure(3);
hold on
plot(wavelength, primaryDist1, 'Color', [1 0 0]);
plot(wavelength, primaryDist2, 'Color', [0 1 0]);
plot(wavelength, primaryDist3, 'Color', [0 0 1]);
xlabel('Wavelength (nm)');
ylabel('Relative Power');
title('Primary (Match) Lights Spectral Distributions')
legend('Primary 1', 'Primary 2', 'Primary 3')
hold off
figure(4);
hold on
plot(wavelength, L_cone, 'Color', [1 0 0]);
plot(wavelength, M_cone, 'Color', [0 1 0]);
plot(wavelength, S_cone, 'Color', [0 0 1]);
xlabel('Wavelength (nm)');
ylabel('Relative Power');
title('Cone Sensitivity');
hold off

%%% Finding the photopic matching system matrix (cone color matching functions) %%%

photopicMatchingSystemMatrix = zeros(3, length(wavelength));

for i = 1:length(wavelength) % Looping through wavelengths
    unitTestLight =  zeros(length(wavelength), 1); 
    unitTestLight(i) = 1; % Set intensity of test light at the current wavelength to 1
    excitations = coneSensitivity' * unitTestLight;
    Q_inv = coneSensitivity' * primaryDist;
    Q = inv(Q_inv);
    column = Q * excitations;
    photopicMatchingSystemMatrix(:, i) = column;
    % We found a linear transformation of the cone spectral sensitivities
end

% Now want to simulate the subject turning knobs %
% The knob settings are within a linear transformation of the cone
% spectral sensitivities
photopicMatchingKnobs = zeros(3, length(wavelength));

% Options for fmincon
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point');

for i = 1:length(wavelength) % Looping through wavelengths

    unitTestLight =  zeros(length(wavelength), 1); 
    unitTestLight(i) = 1; % Set intensity of test light at the current wavelength to 1
    excitations = coneSensitivity' * unitTestLight;

    % Objective function: minimize the difference between the two sides of
    % the equation
    objective = @(knobs) norm(excitations - (coneSensitivity' * (primaryDist * knobs)))^2;

    % Initial guess for the knobs
    initialKnobs = zeros(3, 1);

    % Use fmincon to find knobs for the current wavelength
    % No lower or upper bounds here - add these?
    knobs = fmincon(objective, initialKnobs, [], [], [], [], [], [], [], []);

    photopicMatchingKnobs(:, i) = knobs;
    
end

%%% Plotting photopic matching system matrix in two ways %%%
% They are the same! %
figure(5);
hold on
h1 = plot(wavelength, photopicMatchingKnobs, 'LineWidth', 4, 'LineStyle', '-', 'Color', 'r');
h2 = plot(wavelength, photopicMatchingSystemMatrix, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'b');
legend('Photopic Matching System Matrix Using fmincon', '', '', 'Photopic Matching System Matrix');
xlabel('Wavelength (nm)');
ylabel('Relative Power');
hold off