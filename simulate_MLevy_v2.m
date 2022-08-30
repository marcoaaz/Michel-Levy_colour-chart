%Written by Marco A. Acevedo Zamora, QUT. 26-Aug-2022 
% Following: 
%2013_Sorensen_EJM_Simulating levy chart
%?_2012_Zawischa_Leibniz university Hannover_HOW TO CALCULATE AND RENDER COLOURS IN THIN FILMS

%% Root directory
clear
clc 
close all 

workingDir = 'C:\Users\n10832084\OneDrive - Queensland University of Technology\Desktop\modelling XPL and EBSD\simulating_MLevy';
cd(workingDir)

%% Input

range_fraction = 1000; %chart resolution control

%Sensitivity functions
sensitivityFunctions = readtable('lin2012xyz2e_fine_7sf.csv'); %lambda, R, G, B
sensi_lambda = sensitivityFunctions{:, 1};

sensi_R = sensitivityFunctions{:, 2};
[xData, yData] = prepareCurveData(sensi_lambda, sensi_R);
[fitresult_R, gof_R] = fit(xData, yData, 'linearinterp', 'Normalize', 'on');
sensi_G = sensitivityFunctions{:, 3};
[xData, yData] = prepareCurveData(sensi_lambda, sensi_G);
[fitresult_G, gof_G] = fit(xData, yData, 'linearinterp', 'Normalize', 'on');
sensi_B = sensitivityFunctions{:, 4};
[xData, yData] = prepareCurveData(sensi_lambda, sensi_B);
[fitresult_B, gof_B] = fit(xData, yData, 'linearinterp', 'Normalize', 'on');

%Adobe RGB (Pascale, 2003)
MRGB = [2.04414, -0.5649, -0.3447;
    -0.9693, 1.8760, 0.0416;
    0.0134, -0.1184, 1.0154];

min_lambda = 360;
max_lambda = 830; %nm
min_retardation = 0;
max_retardation = 2500; %nm

thick_var = 30; %microns
% gamma_var = thick_var*birref_var; %retardation
% birref_var = n_gamma - n_alpha;

%Angles:
phi_var = 90*(pi/180); 
tau_var = (45)*(pi/180); 
%phi_var: between polarizer - analyzer =90 
%tau_var: between polarizer - closest 'n' (priviledged?)
%0 degree= XPL-min; 
%45 degree= XPL-max (Michel-Levy diagraM)

%% Simulation

%Spectral transmission
range_retardation = max_retardation - min_retardation;
range_lambda = max_lambda - min_lambda;

gamma_var = min_retardation:range_retardation/range_fraction:max_retardation;
lambda_var = min_lambda:range_lambda/range_fraction:max_lambda; %light wavelength (nm)
n_cols = length(gamma_var);
n_rows = length(lambda_var);

spectral_transmission = zeros(n_rows, n_cols);
k2 = 0 ;
for j = gamma_var
    k2 = k2 + 1;    
    k1 = 0;
    for i = lambda_var
        k1 = k1 + 1;
        
    %Bloss (1999) equation:
    %XPL 
        spectral_transmission(k1, k2) = sin(180*(j/i)*(pi/180))^2; 
    %General
%         spectral_transmission(k1, k2) = (cos(phi_var)^2 - ...
%             (sin(tau_var - phi_var)^2)*...
%             (sin(tau_var)^2)*...
%             (sin(180*(j/i)*(pi/180))^2));     

    end
end

%Human sensitivity
sim_R = fitresult_R(lambda_var); 
sim_G = fitresult_G(lambda_var);
sim_B = fitresult_B(lambda_var);
sim_mat = [sim_R, sim_G, sim_B]';
sim_mat(sim_mat < 0) = 0; %zeroing interpolation

humanVision = rescale(sim_mat*spectral_transmission, 0, 1);
RGB_linear = MRGB*humanVision; %Adobe colour
RGB_linear(RGB_linear > 1) = 1;
RGB_linear(RGB_linear < 0) = 0;
RGB_linear2 = imadjust(RGB_linear, [0, 1], [0, 1], 0.5);
rgb = cat(3, RGB_linear2(1, :), RGB_linear2(2, :), RGB_linear2(3, :));

%% Plots

close all

figure,
plot(lambda_var, sim_R)
hold on
plot(lambda_var, sim_G)
plot(lambda_var, sim_B)
plot(sensi_lambda, sensi_R, 'r+', 'MarkerSize', 1)
plot(sensi_lambda, sensi_G, 'g+', 'MarkerSize', 1)
plot(sensi_lambda, sensi_B, 'b+', 'MarkerSize', 1)
hold off
grid on

figure,
imshow(spectral_transmission, [])

%Michel-Levy chart
chartAspectRatio = 3; %width/height
new_y_max = round(range_fraction/chartAspectRatio);
new_x_max = size(rgb, 2);
x_temp = 0:new_x_max;
rgb_fig = repmat(rgb, [new_y_max, 1]); %image

%Chart customization
max_thickness = 50;
max_retardation = 2500;

real_y = 0:5:max_thickness;
real_x = 0:500:max_retardation; 
y_temp_labels = round(real_y*(new_y_max/max_thickness));
x_temp_labels = round(real_x*(new_x_max/max_retardation));

%prevent Inf (0 <--> 0.0001)
lines_bf = [0.0001, 0.005, 0.01, 0.015, 0.02, ...
    0.025, 0.03, 0.035, 0.04, 0.045, 0.05, ...
    0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2];
n_lines = length(lines_bf);

%text
x_str = strsplit(num2str(real_x));
y_str = strsplit(num2str(real_y));
lines_str = strsplit(num2str(lines_bf));

hFig = figure;
hFig.Position = [50, 100, 1500, 700];

imshow(rgb_fig)
hold on
text_x = zeros(1, n_lines);
text_y = zeros(1, n_lines);
for k = 1:n_lines
    f_temp = x_temp/(lines_bf(k)*(3/0.05));    
    plot(x_temp, f_temp, 'k-');
    
    %text locations
    if lines_bf(k) <= 0.05 %horizontal axis

        [x_data, y_data] = prepareCurveData(f_temp, x_temp);
        [fitresult, ~] = fit(x_data, y_data, 'linearinterp', 'Normalize', 'on');
        text_x(k) = fitresult(new_y_max);
        text_y(k) = new_y_max;

    else %vertical axis

        [x_data, y_data] = prepareCurveData(x_temp, f_temp);
        [fitresult, ~] = fit(x_data, y_data, 'linearinterp', 'Normalize', 'on');
        text_y(k) = fitresult(new_x_max);
        text_x(k) = new_x_max;
    end
    
end
t = text(text_x' + 5, text_y' + 5, lines_str');
% t.Color = 'k';
hold off

%axis settings
ax = gca;
set(ax,'YDir','normal')
grid(ax, "on")
ax.GridColor = 'k';
ax.GridAlpha = 1; 
axis on
xlim([0, new_x_max])
ylim([0, new_y_max])
xticks(x_temp_labels);
yticks(y_temp_labels)
xticklabels(x_str)
yticklabels(y_str)
set(zoom(gcf),'Motion','horizontal','Enable','on');
xlabel('Path difference (nm)')
ylabel('Sample thickness (microns)')
title('Calculated Michel-Levy Colour Chart')
subtitle('Birrefringence')


%%

[fig_rows, fig_cols, n_channels] = size(rgb_fig);
n_pixels = fig_rows*fig_cols;
[row,col] = ind2sub([fig_rows, fig_cols], 1:n_pixels);
x= col';
y= row';

idx_maxBF = (lines_bf == 0.06);
val_top_x = text_x(idx_maxBF); 
val_top_y = text_y(idx_maxBF);
inputMax = val_top_x/val_top_y;
coeff_val = (x./y)*(new_x_max/inputMax);
new_x = coeff_val;

trans_coord = [new_x, y]; %(x, y)
original_coord = [x, y];

r = rgb_fig(:, :, 1); %red   
g = rgb_fig(:, :, 2); %red   
b = rgb_fig(:, :, 3);
r = r(:);
g = g(:);
b = b(:);

%Interpolation
% Warning: Duplicate data points have been detected and removed -
% corresponding values have been averaged.
vec1 = 1; %desired grid
vec2 = 2;
method = 'linear'; %'linear', 'nearest', or 'natural'

F_c = scatteredInterpolant(trans_coord, r, method);
F_c.ExtrapolationMethod = 'linear'; %'nearest', 'linear', or 'none'.
V_r = F_c(original_coord(:, vec1), original_coord(:, vec2)); 

F_c.Values = g; 
V_g = F_c(original_coord(:, vec1), original_coord(:, vec2)); 

F_c.Values = b; 
V_b = F_c(original_coord(:, vec1), original_coord(:, vec2)); 

channels = [V_r, V_g, V_b];        
img_interp = reshape(channels, [fig_rows, fig_cols, n_channels]);

%% Plot

%sampling grid density
figure, 
scatter(trans_coord(:, 1), trans_coord(:, 2), 'r.')
xlim([0, val_top_x])
ylim([0, val_top_y])
grid on
axis equal

% Raith–Sørensen chart

%Chart customization
max_bifferingence = 0.06;
real_x = 0:0.002:max_bifferingence; 
% real_x_lab = 0:0.01:max_bifferingence; 
x_temp_labels = round(real_x*(new_x_max/max_bifferingence));

%prevent Inf (0 <--> 0.0001)
lines_bf = [0.0001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06];
n_lines = length(lines_bf);

%text
x_str = strsplit(num2str(real_x));

close all

hFig2 = figure;
hFig2.Position = [50, 100, 1500, 700];

imshow(img_interp)
ax = gca;
set(ax, 'YDir', 'normal')
grid(ax, "on")
ax.GridColor = 'k';
ax.GridAlpha = 0.3; 

%axis settings
ax = gca;
set(ax,'YDir','normal')
grid(ax, "on")
ax.GridColor = 'k';
ax.GridAlpha = 1; 
axis on
xlim([0, new_x_max])
ylim([0, new_y_max])
xticks(x_temp_labels);
yticks(y_temp_labels)
xticklabels(x_str)
yticklabels(y_str)
set(zoom(gcf), 'Motion', 'horizontal', 'Enable','on');
xlabel('Birrefringence')
ylabel('Sample thickness (microns)')
title('Raith–Sørensen chart')


