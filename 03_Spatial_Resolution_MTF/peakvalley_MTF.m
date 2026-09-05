%% MTF Calculator with Mean Calculation
% Υπολογίζει το MTF παίρνοντας τον Μέσο Όρο (Mean) πολλών κορυφών/κοιλάδων.
clc; clear; close all;

% =========================================================================
% 1. ΕΙΣΑΓΩΓΗ ΜΕΤΡΗΣΕΩΝ 
% =========================================================================
% Οδηγίες:
% Για κάθε συχνότητα, φτιάχνουμε μια "δομή" (structure).
% data(i).freq    -> Η συχνότητα (lp/cm)
% data(i).peaks   -> Ολες οι τιμές των ΚΟΡΥΦΩΝ που μετρησες (μες στις αγκύλες)
% data(i).valleys -> Ολες οι τιμές των ΚΟΙΛΑΔΩΝ που μετρησες

% --- ΜΕΤΡΗΣΗ 1 (π.χ. 1 lp/cm) ---
i = 1;
data(i).freq = 1; 
data(i).peaks = [2775, 2749];   % Βάλε εδώ όσες κορυφές βρήκες
data(i).valleys = [-60];         % Βάλε εδώ όσες κοιλάδες βρήκες

% --- ΜΕΤΡΗΣΗ 2 (π.χ. 3 lp/cm) ---
i = i + 1;
data(i).freq = 2;
data(i).peaks = [2593, 2519, 2505];   % Παράδειγμα: 3 κορυφές
data(i).valleys = [38, 109];    % Παράδειγμα: 3 κοιλάδες

% --- ΜΕΤΡΗΣΗ 3 (π.χ. 5 lp/cm) ---
i = i + 1;
data(i).freq = 3;
data(i).peaks = [2195, 2081, 1921, 1870];         
data(i).valleys = [243, 371, 435];

% --- ΜΕΤΡΗΣΗ 4 (π.χ. 7 lp/cm) ---
i = i + 1;
data(i).freq = 4;
data(i).peaks = [1527, 1515, 1543, 1563];         
data(i).valleys = [812, 704, 629];

% --- ΜΕΤΡΗΣΗ 5 (π.χ. 7 lp/cm) ---
i = i + 1;
data(i).freq = 5;
data(i).peaks = [1100, 1180, 1107, 1206];         
data(i).valleys = [891, 849, 866];

% --- ΜΕΤΡΗΣΗ 6 (π.χ. 7 lp/cm) ---
i = i + 1;
data(i).freq = 6;
data(i).peaks = [1021, 1041, 997, 1164, 1028];         
data(i).valleys = [860, 905, 880, 923];

% --- ΜΕΤΡΗΣΗ 7 (π.χ. 7 lp/cm) ---
i = i + 1;
data(i).freq = 7;
data(i).peaks = [984, 1023, 1007, 992];         
data(i).valleys = [973, 958];

% --- ΜΕΤΡΗΣΗ 8 (π.χ. 7 lp/cm) ---
i = i + 1;
data(i).freq = 8;
data(i).peaks = [1061, 1088, 1038];         
data(i).valleys = [1012];

% --- ΜΕΤΡΗΣΗ 9 (π.χ. 7 lp/cm) ---
i = i + 1;
data(i).freq = 9;
data(i).peaks = [1126, 1119];         
data(i).valleys = [1103];
% (Αν θες κι άλλες, κάνε copy-paste το μπλοκ "i = i + 1" από πάνω)

% =========================================================================
% 2. ΥΠΟΛΟΓΙΣΜΟΙ (Μην πειράξεις τίποτα από εδώ και κάτω)
% =========================================================================

num_measurements = length(data);
results_freq = zeros(num_measurements, 1);
results_mod = zeros(num_measurements, 1);
mean_peaks = zeros(num_measurements, 1);
mean_valleys = zeros(num_measurements, 1);

fprintf('\n--- ΑΝΑΛΥΤΙΚΑ ΑΠΟΤΕΛΕΣΜΑΤΑ ---\n');
fprintf('%-10s %-12s %-12s %-12s\n', 'Freq', 'Mean Peak', 'Mean Valley', 'Modulation');

for k = 1:num_measurements
    % 1. Υπολογισμός Μέσων Όρων
    m_peak = mean(data(k).peaks);
    m_valley = mean(data(k).valleys);
    
    % 2. Υπολογισμός Modulation: (Max-Min)/(Max+Min)
    % Προσθέτουμε +1000 αν οι τιμές είναι σε HU (ώστε να αποφύγουμε το 0)
    % Αλλά αν οι τιμές σου είναι ήδη θετικές (π.χ. 2775), ο τύπος δουλεύει απευθείας.
    mod_val = (m_peak - m_valley) / (m_peak + m_valley);
    
    % Αποθήκευση
    results_freq(k) = data(k).freq;
    results_mod(k) = mod_val;
    mean_peaks(k) = m_peak;
    mean_valleys(k) = m_valley;
    
    fprintf('%-10.1f %-12.1f %-12.1f %-12.4f\n', data(k).freq, m_peak, m_valley, mod_val);
end

% 3. Κανονικοποίηση (MTF)
% Θεωρούμε την χαμηλότερη συχνότητα ως το 100% (Reference)
[sorted_freq, sort_idx] = sort(results_freq);
sorted_mod = results_mod(sort_idx);

mtf_curve = sorted_mod / sorted_mod(1);

%% =========================================================================
% 3. ΓΡΑΦΗΜΑ & ΣΥΜΠΕΡΑΣΜΑΤΑ
% =========================================================================

% Πίνακας Τελικός
T = table(sorted_freq, sorted_mod, mtf_curve, ...
    'VariableNames', {'Frequency_lp_cm', 'Raw_Modulation', 'MTF_Normalized'});
disp(' ');
disp('--- ΤΕΛΙΚΟΣ ΠΙΝΑΚΑΣ CTF ---');
disp(T);

% Γράφημα
figure('Name', 'CTF Curve Analysis', 'Color', 'w');
plot(sorted_freq, mtf_curve, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
grid on;
xlabel('Spatial Frequency (lp/cm)', 'FontSize', 12);
ylabel('Modulation Transfer Function (MTF)', 'FontSize', 12);
title('MTF Curve (Average Peak/Valley Method)', 'FontSize', 14);
ylim([0 1.1]);

% Interpolation για 50% και 10%
try
    mtf50 = interp1(mtf_curve, sorted_freq, 0.5, 'pchip');
    mtf10 = interp1(mtf_curve, sorted_freq, 0.1, 'pchip');
    
    yline(0.5, '--r', ['CTF 50%: ' num2str(mtf50, '%.2f') ' lp/cm']);
    yline(0.1, '--m', ['CTF 10%: ' num2str(mtf10, '%.2f') ' lp/cm']);
    
    fprintf('\nCTF 50%% (50%% Contrast): %.2f lp/cm\n', mtf50);
    fprintf('CTF 10%% (Limiting Resolution): %.2f lp/cm\n', mtf10);
catch
    disp('Δεν υπάρχουν αρκετά σημεία για ακριβή υπολογισμό MTF50/10.');
end