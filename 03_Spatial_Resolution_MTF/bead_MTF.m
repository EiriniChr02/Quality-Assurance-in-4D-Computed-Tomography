%% MTF Calculator using Gaussian Fit (The Robust Method)
% Αυτή η μέθοδος είναι ιδανική όταν έχουμε λίγα σημεία (Bead/Wire).
% Προσαρμόζει μια καμπύλη Gauss στα δεδομένα για τέλειο MTF χωρίς θόρυβο.

clc; clear; close all;

% =========================================================================
% 1. ΕΙΣΑΓΩΓΗ ΔΕΔΟΜΕΝΩΝ
% =========================================================================

prompt = {'Κάνε Paste τα δεδομένα από το ImageJ:', 'Pixel Size (mm) (π.χ. 0.5):'};
dlgtitle = 'Εισαγωγή Δεδομένων';
dims = [20 80; 1 50];
definput = {'', '0.8'}; 
answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer), return; end

% Επεξεργασία εισόδου
raw_text = answer{1};
pixel_size_mm = str2double(answer{2});
data_matrix = str2num(raw_text);

if isempty(data_matrix)
    errordlg('Δεν βρέθηκαν αριθμοί.', 'Σφάλμα'); return;
end

% Επιλογή στήλης (αν έχεις 2 στήλες, παίρνει τη 2η)
[~, cols] = size(data_matrix);
if cols >= 2
    y_raw = data_matrix(:, 2); 
else
    y_raw = data_matrix(:, 1);
end
y_raw = y_raw(:)'; % Μετατροπή σε γραμμή
x_raw = 1:length(y_raw); % Ο άξονας των pixels

% =========================================================================
% 2. ΠΡΟΕΤΟΙΜΑΣΙΑ & GAUSSIAN FITTING
% =========================================================================

% Α. Σωστή Αφαίρεση Υποβάθρου (Background Subtraction)
% Παίρνουμε τον μέσο όρο των 2 πρώτων και 2 τελευταίων σημείων ως background
bkg = mean([y_raw(1:2), y_raw(end-1:end)]);
y_clean = y_raw - bkg;
y_clean(y_clean < 0) = 0; % Απαλοιφή αρνητικών τιμών

% Β. Προσαρμογή Καμπύλης Gauss (Gaussian Fit)
% Μοντέλο: y = A * exp( -((x-mu)^2) / (2*sigma^2) )
% Χρησιμοποιούμε fminsearch για να βρούμε τα A, mu, sigma που ταιριάζουν τέλεια

% Αρχικές προβλέψεις (Guesses)
[max_val, max_idx] = max(y_clean);
start_params = [max_val, max_idx, 1.0]; % A=max, mu=index, sigma=1 pixel

% Συνάρτηση σφάλματος (για ελαχιστοποίηση)
gauss_fun = @(p, x) p(1) * exp( -((x - p(2)).^2) / (2*p(3)^2) );
cost_fun = @(p) sum((y_clean - gauss_fun(p, x_raw)).^2);

% Εκτέλεση Fitting
params = fminsearch(cost_fun, start_params);

A_fit = params(1);
mu_fit = params(2);
sigma_fit = abs(params(3)); % Το sigma πρέπει να είναι θετικό

% =========================================================================
% 3. ΔΗΜΙΟΥΡΓΙΑ ΥΨΗΛΗΣ ΑΝΑΛΥΣΗΣ & FFT
% =========================================================================

% Δημιουργούμε μια "Τέλεια" καμπύλη με πολύ υψηλή ανάλυση (Oversampling)
% Αυτό αντικαθιστά το απλό Zero Padding και είναι πιο φυσικό.
step_new = 0.1; % Υπο-pixel ανάλυση (10 φορές πυκνότερη)
x_high_res = 1 : step_new : length(y_raw);
lsf_fitted = gauss_fun(params, x_high_res);

% FFT στην προσαρμοσμένη καμπύλη
N_fft = 2048; % Μεγάλο νούμερο για ομαλότητα
fft_res = abs(fft(lsf_fitted, N_fft));

% MTF Curve
num_points = floor(N_fft/2);
mtf_curve = fft_res(1:num_points);
mtf_curve = mtf_curve / mtf_curve(1); % Κανονικοποίηση στο 1.0

% Άξονας Συχνοτήτων
% Προσοχή: Το sampling rate τώρα είναι 1 / (pixel_size * step_new)
fs = 1 / (pixel_size_mm * step_new); % samples/mm
fs_lp_cm = fs * 10; % lp/cm
df = fs_lp_cm / N_fft;
freq_axis = (0:num_points-1) * df;

% =========================================================================
% 4. ΑΠΟΤΕΛΕΣΜΑΤΑ & ΓΡΑΦΗΜΑΤΑ
% =========================================================================

figure('Name', 'Gaussian Fit MTF', 'Color', 'w');

% Plot 1: Τα δεδομένα σου vs Η Προσαρμογή (Fitting)
subplot(2,1,1);
plot(x_raw, y_clean, 'bo', 'MarkerSize', 6, 'LineWidth', 1.5); hold on;
plot(x_high_res, lsf_fitted, 'r-', 'LineWidth', 2);
legend('Raw Data (ImageJ)', 'Gaussian Fit');
title(['Gaussian Fit on LSF (Sigma = ' num2str(sigma_fit, '%.2f') ' pixels)']);
xlabel('Pixel Position'); ylabel('Intensity (Background Corrected)');
grid on;

% Plot 2: Η Τελική MTF
subplot(2,1,2);
plot(freq_axis, mtf_curve, 'b-', 'LineWidth', 2);
title('MTF Curve (from Gaussian Fit)');
xlabel('Spatial Frequency (lp/cm)'); ylabel('MTF');
xlim([0 20]); ylim([0 1.1]);
grid on;

% Υπολογισμός Τιμών
try
    mtf50 = interp1(mtf_curve, freq_axis, 0.5, 'pchip');
    mtf10 = interp1(mtf_curve, freq_axis, 0.1, 'pchip');
    
    yline(0.5, '--r', ['MTF 50%: ' num2str(mtf50, '%.2f')]);
    yline(0.1, '--m', ['MTF 10%: ' num2str(mtf10, '%.2f')]);
    
    fprintf('\n--- ΑΠΟΤΕΛΕΣΜΑΤΑ (Gaussian Method) ---\n');
    fprintf('MTF 50%%: %.2f lp/cm\n', mtf50);
    fprintf('MTF 10%%: %.2f lp/cm\n', mtf10);
    fprintf('Sigma (Width): %.2f pixels\n', sigma_fit);
catch
    disp('Σφάλμα στον υπολογισμό των τιμών.');
end