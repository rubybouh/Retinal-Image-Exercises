clear all
close all

% Calculating the retinal image using Fourier transforms
nPixels = 101;

% Choose parameters and range of x values
mu = 0;
sigma = 1;
x = linspace(-5, 5, nPixels);

% Compute the values for a symmetric pointspread
sym_pointspread = normpdf(x, mu, sigma);

% POINTSPREAD SELECTION %
pointspread = sym_pointspread;
pointspread = pointspread / sum(pointspread);

% Plot the pointspread
figure(1)
plot(x, pointspread);
xlabel('Retinal Position Unit');
ylabel('Intensity');
title('Pointspread Function');

% Initialize system matrix
system_matrix = zeros(nPixels, nPixels);

% Populate central column with pointspread function values
central_column_index = ceil(nPixels/2);
system_matrix(:, central_column_index) = pointspread;

% Filling in the system matrix

for col_index = 1:nPixels
    if col_index <= central_column_index 
        col_shift = central_column_index - col_index; 
        shifted_col = circshift(pointspread, -col_shift);
    else
        col_shift = col_index - central_column_index;
        shifted_col = circshift(pointspread, col_shift);
    end
    system_matrix(:, col_index) = shifted_col;
end

% To visualize the system matrix
figure; imagesc(system_matrix); axis('square');

% Harmonic Input
A = 2; % Amplitude
f = 0.3; % Frequency
alpha = pi/4; % phase shift

display_image_sine = (A * sin(2 * pi * f * x + alpha))';
display_image_sine = zeros(size(display_image_sine));
display_image_sine(central_column_index) = 1;

figure(2) 
sgtitle('Shift-Invariant Optical System with Harmonic Input');

subplot(2, 1, 1);
plot(x, display_image_sine);
xlabel('Spatial Position Unit');
ylabel('Intensity');
title('Display Image Intensity (Input)');

% Finding retinal image through matrix multiplication and using Fourier
% transforms

retinal_image = system_matrix * display_image_sine;

% Finding and sorting eigenvectors and eigenvalues
[eigenvectors, eigenvalues_matrix] = eig(system_matrix);
% Returns diagonal matrix D of eigenvalues and matrix V whose columns are the corresponding right eigenvectors, so that A*V = V*D.

% Extracting eigenvalues from diagonal matrix
eigenvalues = diag(eigenvalues_matrix);

% Sort the eigenvalues in ascending order and get the sorting indices
[~, sorting_indices] = sort(eigenvalues);

% Reorder the eigenvectors according to the sorted eigenvalues
sorted_eigenvalues = eigenvalues(sorting_indices);
sorted_eigenvectors = eigenvectors(:, sorting_indices);

% Flipping evectors and evalues so evalues are in descending order
D = diag(sorted_eigenvalues(end:-1:1));
desc_eigenvectors = sorted_eigenvectors(:,end:-1:1);

% Verify the reconstruction
reconstructed_system_matrix = desc_eigenvectors * D * inv(desc_eigenvectors);
dist = norm(reconstructed_system_matrix - system_matrix);

% Retinal image formula without shifts
% retinal_image_fft = ifft(diag(real(fft(pointspread)))*fft(display_image_sine));
% Correct formula with shifts
retinal_image_fft = ifftshift(ifft(diag(fft(fftshift(pointspread)))*fft(fftshift(display_image_sine))));

% distances = abs(retinal_image - retinal_image_fft);

% Alternate method: multiply E inverse, SM, and E to get D (D_direct = D)
% D_direct = eigenvectors * system_matrix * inv(eigenvectors);
% retinal_image_fft = inv(eigenvectors) * (D_direct * (eigenvectors * display_image_sine));

% Plotting retinal image with harmonic input
subplot(2, 1, 2);
hold on
plot(x, retinal_image);
plot(x, retinal_image_fft); 
xlabel('Retinal Position Unit');
ylabel('Intensity');
title('Intensity of Retinal Image (Output)');
legend('Matrix Mult.', 'Fourier Transforms')
hold off
% If they don't match up exactly it is due to numerical precision errors
% and/or undersampling

% FINDING RETINAL IMAGE USING MATLAB'S IFFT MATRIX %

% Creating the E_ifft matrix
matrix = diag(ones(nPixels, 1));
E_ifft = zeros(nPixels, nPixels);
E_ifft_reordered = zeros(nPixels, nPixels);

for i = (1:nPixels)
    col = ifftshift(ifft(matrix(:,i)));
    E_ifft(:, i) = col;
end

% Checking that the columns of E_ifft (complex exponentials) correspond to
% linear combinations of columns in E_M

% This did not work, most likely due to numerical precision errors
% quotient1 = zeros(nPixels, central_column_index);
% quotient2 = zeros(nPixels, central_column_index);

% frequency = i;
% eCheck1 = E_M(:,1+2*frequency-1:1+2*frequency)*(E_M(:,1+2*frequency-1:1+2*frequency)\E_ifft(:,1+frequency));
% % We should have eCheck1 = E_ifft(:,1+frequency)
% eCheck2 = E_M(:,1+2*frequency-1:1+2*frequency)*(E_M(:,1+2*frequency-1:1+2*frequency)\E_ifft(:,nPixels-frequency+1));
% % We should have eCheck2 = E_ifft(:,nPixels-frequency+1)
% quotient1(:, i) = eCheck1 .\ E_ifft(:,1+frequency);
% quotient2(:, i) = eCheck2 .\ E_ifft(:,nPixels-frequency+1);

frequencies = 1:(nPixels - 1)/2; 
E_M = desc_eigenvectors;
list1 = [];
list2 = [];

for i = (1:length(frequencies))
    frequency = i;
    % Get the regression weights on the whole eigenvector matrix E_M.
    % Only two should be very big for each column of the E_ifft matrix.
    weights1 = E_M\E_ifft(:,1+frequency);
    weights2 = E_M\E_ifft(:,nPixels-frequency+1);
    % weights1 and weights2 should be the same

    absWeights1 = abs(weights1);
    absWeights2 = abs(weights2);

    % Finding largest and second largest weights
    max_weight1 = max(absWeights1); 
    max_index_weight1 = find(absWeights1 == max_weight1);
    % Removing the max value so that we can find the second largest value
    absWeights1(max_index_weight1) = [];
    second_max_weight1 = max(absWeights1);
    second_max_index_weight1 = find(absWeights1 == second_max_weight1);
    list1{end + 1} = {max(max_index_weight1); max(second_max_index_weight1)};
    % Put max() to account for the weight being the same at two indices
    % (choose the larger one)

    max_weight2 = max(absWeights2); 
    max_index_weight2 = find(absWeights2 == max_weight2);
    % Removing the max value so that we can find the second largest value
    absWeights2(max_index_weight2) = [];
    second_max_weight2 = max(absWeights2);
    second_max_index_weight2 = find(absWeights2 == second_max_weight2);
    list2{end + 1} = {max(max_index_weight2); max(second_max_index_weight2)};
    % Put max() to account for the weight being the same at two indices
    % (choose the larger one)

    % Then reorder the E_ifft matrix
    % Want to move the E_ifft column to the index of the weight
    % Using the weights1 list
    E_ifft_reordered(:,1+frequency) = E_ifft(:, list2{frequency}{1});
    E_ifft_reordered(:,nPixels-frequency+1) = E_ifft(:, (list2{frequency}{1} + 1));
end

% absWeights1 = sort(abs(weights1),'descend');
% absWeights2 = sort(abs(weights2),'descend');

% Using list 2 to reorder the MATLAB ifft matrix - should correspond to
% increasing spatial frequency with increasing index.

D_direct = inv(E_ifft) * system_matrix * E_ifft;
retinal_image_fft_2 = E_ifft * (D_direct * (inv(E_ifft) * display_image_sine));

figure(3) 
sgtitle('Shift-Invariant Optical System with Harmonic Input, using MATLABs IFFT Matrix');

subplot(2, 1, 1);
plot(x, display_image_sine);
xlabel('Spatial Position Unit');
ylabel('Intensity');
title('Display Image Intensity (Input)');

% Plotting retinal image with harmonic input, using E_ifft
subplot(2, 1, 2);
hold on
plot(x, retinal_image);
plot(x, retinal_image_fft_2); 
xlabel('Retinal Position Unit');
ylabel('Intensity');
title('Intensity of Retinal Image (Output)');
legend('Matrix Mult.', 'Fourier Transforms')
hold off

% FREQUENCY PLOT %

figure(4);

fs = 1/0.1;
% 0.1 seconds between each sample (sampling period), so the signal is sampled 10 times per
% second
N = length(display_image_sine);
Freq = (0:N-1)*(fs/N);
% Scaling the indices to frequency values

Y = E_ifft_reordered * (pointspread)';
% Fourier transform of the pointspread

% plot(Freq(1:round(N/2)), Y(1:round(N/2)));
plot(Freq, Y)