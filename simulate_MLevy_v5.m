%Written by Marco A. Acevedo Zamora, QUT. 26-Aug-2022 
% Following: 
%2013_Sorensen_EJM_Simulating levy chart
%?_Zawischa_Leibniz university Hannover_HOW TO CALCULATE AND RENDER COLOURS IN THIN FILMS

%Root directory
clear
clc 
close all 

workingDir = 'E:\Alienware_March 22\current work\00-new code May_22\simulating_MLevy';
cd(workingDir)

scriptDir = 'E:\Alienware_March 22\current work\00-new code May_22\';
addpath(fullfile(scriptDir, '\spectral_color_1'))

%% Input

%Lightness sensitivity functions of the eye (colour matching) :
%Source: http://cvrl.ioo.ucl.ac.uk/index.htm
%Download site: New CIE XYZ functions transformed from the CIE (2006) LMS functions
%Download file: 2-deg XYZ CMFs transformed from the CIE (2006) 2-deg LMS cone fundamentals
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

%Adobe RGB colour model (Pascale, 2003): transformation to digital device display
MRGB = [2.04414, -0.5649, -0.3447;
    -0.9693, 1.8760, 0.0416;
    0.0134, -0.1184, 1.0154];

%ranges for estimation (according to charts)
min_lambda = 360;
max_lambda = 830; %nm
min_retardation = 0;
max_retardation = 2500; %nm

%General equation parameters (dont worry about this)
thick_var = 30; %microns
% birref_var = n_gamma - n_alpha;
% gamma_var = thick_var*birref_var; %retardation

%Angles:
phi_var = 90*(pi/180); 
tau_var = (45)*(pi/180); 
%phi_var: between polarizer - analyzer =90 
%tau_var: between polarizer - closest 'n' (priviledged?)
%0 degree= XPL-min; 
%45 degree= XPL-max (only max birrefringence plot in Michel-Levy diagram)

%User input:
chartAspectRatio = 3; %width/height
range_fraction = 1000; %chart resolution control
chart_x_min_real = 0;
chart_x_max_real = max_retardation;
chart_y_min_real = 30; %microns
chart_y_max_real = 35;
max_thickness = 50; %published in books

%% Simulation

%Spectral transmission
range_retardation = max_retardation - min_retardation;
range_lambda = max_lambda - min_lambda;

gamma_var = min_retardation:range_retardation/range_fraction:max_retardation;
lambda_var = min_lambda:range_lambda/range_fraction:max_lambda; %light wavelength (nm)
n_cols = length(gamma_var);
n_rows = length(lambda_var);

spectral_transmission = zeros(n_rows, n_cols); %colour stimulus
k2 = 0 ;
for j = gamma_var %retardation
    k2 = k2 + 1;    
    k1 = 0;
    for i = lambda_var %wavelength
        k1 = k1 + 1;
        
    %Bloss (1999) equation: The variation of the intensity of the light
    %being transmitted by the analyzer      
    % XPL
        spectral_transmission(k1, k2) = sin(180*(j/i)*(pi/180))^2; 

    %General
%         spectral_transmission(k1, k2) = (cos(phi_var)^2 - ...
%             (sin(tau_var - phi_var)^2)*...
%             (sin(tau_var)^2)*...
%             (sin(180*(j/i)*(pi/180))^2));     

    end
end

%Human sensitivity correction
sim_R = fitresult_R(lambda_var); 
sim_G = fitresult_G(lambda_var);
sim_B = fitresult_B(lambda_var);
sim_mat = [sim_R, sim_G, sim_B]'; %(3 x Lambdas)
sim_mat(sim_mat < 0) = 0; %zeroing interpolation
humanVision = rescale(sim_mat*spectral_transmission, 0, 1); %(3 x Retardations)

%Display device correction
RGB_linear = MRGB*humanVision; %Adobe colour
RGB_linear(RGB_linear > 1) = 1;
RGB_linear(RGB_linear < 0) = 0;
RGB_linear2 = imadjust(RGB_linear, [0, 1], [0, 1], 0.5); %gamma correction (LCD monitor)
rgb = cat(3, RGB_linear2(1, :), RGB_linear2(2, :), RGB_linear2(3, :)); %1 strip of pixels

%% Plot 0 (optional): Chromacity diagram (test)
%Following: 2012_Zawischa_How to calculate and render colours in thin films
%oil bubble double reflection interference, D65 sky white source

testVision = sim_mat*spectral_transmission; %XYZ to be converted in tristimulus values (RGB)
x_colour = testVision(1, :)./sum(testVision, 1);
y_colour = testVision(2, :)./sum(testVision, 1);

xy_white = [0.3127, 0.3290];
xy_R = [.64, 0.33];
xy_G = [.3, 0.6];
xy_B = [.15, 0.06];
display_gamma = 1/2.4; %CRT monitor= 2.4

clear sol u1 u2 u3 Rx Gx Bx Ry Gy By Rz Gz Bz
syms u1 u2 u3 Rx Gx Bx Ry Gy By Rz Gz Bz

XYZ_white = [xy_white(1)/xy_white(2), 1, (1-xy_white(1)-xy_white(2))/xy_white(2)];

equ1 = xy_white(1)*u1 == XYZ_white(1);
equ2 = xy_white(2)*u2 == XYZ_white(2);
equ3 = (1-xy_white(1)-xy_white(2))*u3 == XYZ_white(3);
sol1 = solve([equ1, equ2, equ3], [u1, u2, u3]);

clear u1 u2 u3 

u1 = double(sol1.u1);
u2 = double(sol1.u2);
u3 = double(sol1.u3);

XYZ_R = u1*[xy_R(1), xy_R(2), 1-xy_R(1)-xy_R(2)];
XYZ_G = u2*[xy_G(1), xy_G(2), 1-xy_G(1)-xy_G(2)];
XYZ_B = u3*[xy_B(1), xy_B(2), 1-xy_B(1)-xy_B(2)];
XYZ_system = [XYZ_R; XYZ_G; XYZ_B].';
X_vect = [Rx; Gx; Bx];
Y_vect = [Ry; Gy; By];
Z_vect = [Rz; Gz; Bz];

eq1 = XYZ_system(1, :)*X_vect == 1;
eq2 = XYZ_system(2, :)*X_vect == 0;
eq3 = XYZ_system(3, :)*X_vect == 0;
eq4 = XYZ_system(1, :)*Y_vect == 0;
eq5 = XYZ_system(2, :)*Y_vect == 1;
eq6 = XYZ_system(3, :)*Y_vect == 0;
eq7 = XYZ_system(1, :)*Z_vect == 0;
eq8 = XYZ_system(2, :)*Z_vect == 0;
eq9 = XYZ_system(3, :)*Z_vect == 1;

sol = solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9], ...
    [Rx, Gx, Bx, Ry, Gy, By, Rz, Gz, Bz]);

%range for simulation
step = 5;
X_sim = 0:step:255; %edit minimum to 1
Y_sim = 0:step:255;
Z_sim = 0:step:255;
[X_sim_mat, Y_sim_mat, Z_sim_mat] = meshgrid(X_sim, Y_sim, Z_sim);
[n_rows, n_cols, n_channels] = size(X_sim_mat);
n_spots = n_rows*n_cols*n_channels;
X_sim_vec = X_sim_mat(:);
Y_sim_vec = Y_sim_mat(:);
Z_sim_vec = Z_sim_mat(:);

testVision_sim = [X_sim_vec, Y_sim_vec, Z_sim_vec].';
x_colour_sim = testVision_sim(1, :)./sum(testVision_sim, 1);
y_colour_sim = testVision_sim(2, :)./sum(testVision_sim, 1);

Y_R = u1*xy_R(2);
Y_G = u2*xy_G(2);
Y_B = u3*xy_B(2);
factor = [Y_R; Y_G; Y_B];

%vectors
R = double(sol.Rx)*X_sim_vec + double(sol.Ry)*Y_sim_vec + double(sol.Rz)*Z_sim_vec;
G = double(sol.Gx)*X_sim_vec + double(sol.Gy)*Y_sim_vec + double(sol.Gz)*Z_sim_vec;
B = double(sol.Bx)*X_sim_vec + double(sol.By)*Y_sim_vec + double(sol.Bz)*Z_sim_vec;

max_value = 0;
for i = 1:n_spots
    tristimulus_temp = [R(i), G(i), B(i)];
    max_temp = max(tristimulus_temp);
    max_value = max(max_temp, max_value);
end
tristimulus_vec = [R, G, B]./max_value;

%checking range
figure

histogram(tristimulus_vec(:, 2), 'FaceColor', 'g', 'FaceAlpha', 1)
hold on
histogram(tristimulus_vec(:, 1), 'FaceColor', 'r', 'FaceAlpha', 1)
histogram(tristimulus_vec(:, 3), 'FaceColor', 'b', 'FaceAlpha', 1)
hold off
grid on
legend({'g', 'r', 'b'})

%Correcting tristimulus
tristimulus_corr = zeros(n_spots, 3);
brightness = zeros(n_spots, 1);
luminance = zeros(n_spots, 1);
for i = 1:n_spots     %n_spots    
    tristimulus_temp = tristimulus_vec(i, :);
    idx_negative = (tristimulus_temp < 0);
    
    factor_vec = factor.*(~idx_negative).';
    brightness_orig = tristimulus_temp*factor;
    brightness_temp = tristimulus_temp*(factor_vec); %brightness
%     brightness_temp = [R, G, B]*(factor); %brightness

    correction_vec = brightness_temp*[
        1/(brightness_temp - factor(1)*tristimulus_temp(1)); 
        1/(brightness_temp - factor(2)*tristimulus_temp(2)); 
        1/(brightness_temp - factor(3)*tristimulus_temp(3)); 
        ];
    
    array_temp = correction_vec(~idx_negative);
    correction_factor = prod(array_temp, 'all' );

    %correct luminance
    brightness_temp2 = tristimulus_temp*(factor_vec.*correction_factor); 
    tristimulus_temp2 = correction_factor*(tristimulus_temp.*(~idx_negative)); 

    %assigning
    tristimulus_corr(i, :) = tristimulus_temp2;
    brightness(i) = brightness_orig;
    luminance(i) = brightness_temp2;

%     disp(i)
end

%sRGB conversion
idx_neg_points = (luminance <= 0);
% idx_outlier = (luminance == max(luminance));
idx_nan_points = isnan(luminance);

%filter luminance
% idx_out = idx_neg_points | idx_nan_points | idx_outlier;
idx_out = idx_neg_points | idx_nan_points;
n_posi = sum(~idx_out);
x_posi = x_colour_sim(~idx_out); %
x_nega = x_colour_sim(idx_out);
y_posi = y_colour_sim(~idx_out); %
y_nega = y_colour_sim(idx_out);
tri_posi = tristimulus_corr(~idx_out, :); %linear
 
tri_sRGB = zeros(n_posi, 3); %non-linear birghtness perception
for k = 1:n_posi
    array_temp = tri_posi(k, :);
    for m = 1:3
        if array_temp(m) <= 0.00304
            tri_sRGB(k, m) = array_temp(m)*12.92;
        else
            tri_sRGB(k, m) = 1.055*(array_temp(m)^display_gamma)-0.055;
        end        
    end
end
% tri_sRGB2 = uint8(rescale(tri_sRGB, 0, 255, 'InputMax', 6)); %manual edit
tri_sRGB2 = 255*tri_sRGB;

%triangle filter
idx_0 = (tri_sRGB2(:, 1) == 0) | (tri_sRGB2(:, 2) == 0) | (tri_sRGB2(:, 3) == 0);
idx_1 = (tri_sRGB2(:, 1) > 255) & (tri_sRGB2(:, 2) > 255) & (tri_sRGB2(:, 3) > 255);
idx_2 = (tri_sRGB2(:, 1) == 0) & ((tri_sRGB2(:, 2) > 255) | (tri_sRGB2(:, 3) > 255));
idx_3 = (tri_sRGB2(:, 2) == 0) & ((tri_sRGB2(:, 1) > 255) | (tri_sRGB2(:, 3) > 255));
idx_4 = (tri_sRGB2(:, 3) == 0) & ((tri_sRGB2(:, 1) > 255) | (tri_sRGB2(:, 2) > 255));
idx_oversaturated = idx_0 | idx_1 | idx_2 | idx_3 | idx_4;

tri_sRGB3 = tri_sRGB2(~idx_oversaturated, :); %not converted to uint8
x_posi2 = x_posi(~idx_oversaturated);
y_posi2 = y_posi(~idx_oversaturated);
%check for filtering
trans_test = [x_posi2', y_posi2', tri_sRGB3];

r = double(tri_sRGB3(:, 1));
g = double(tri_sRGB3(:, 2));
b = double(tri_sRGB3(:, 3));

%generating triangle image
new_x_max = 0.8;
new_y_max = 0.9;
new_x = 0:1/range_fraction:new_x_max;
new_y = 0:1/range_fraction:new_y_max; %light wavelength (nm)
n_cols1 = length(new_x);
n_rows1 = length(new_y);
[new_x_m, new_y_m] = meshgrid(new_x, new_y);
new_x_v = new_x_m(:);
new_y_v = new_y_m(:);

original_coord = [new_x_v, new_y_v];
trans_coord = [x_posi2', y_posi2'];
[trans_coord_u, i_unique, ~] = unique(trans_coord, 'rows');
r_u = r(i_unique);
g_u = g(i_unique);
b_u = b(i_unique);

%Interpolation
% Warning: Duplicate data points have been detected and removed -
% corresponding values have been averaged. 
% The number of data point locations should equal the number of data point values.

vec1 = 1; %desired grid
vec2 = 2;
method = 'linear'; %'linear', 'nearest', or 'natural'

F_c = scatteredInterpolant(trans_coord_u, r_u, method);
F_c.ExtrapolationMethod = 'none'; %'nearest', 'linear', or 'none'.
V_r = F_c(original_coord(:, vec1), original_coord(:, vec2)); 

F_c.Values = g_u; 
V_g = F_c(original_coord(:, vec1), original_coord(:, vec2)); 

F_c.Values = b_u; 
V_b = F_c(original_coord(:, vec1), original_coord(:, vec2)); 

channels = [V_r, V_g, V_b];  

% img_interp = uint8(reshape(channels, [n_rows1, n_cols1, 3])); %paper method
img_interp_rs = uint8(rescale(reshape(channels, [n_rows1, n_cols1, 3]), 0, 255, 'InputMax', 500));
%patch for visualization convenience
r_temp = img_interp_rs(:, :, 1);
g_temp = img_interp_rs(:, :, 2);
b_temp = img_interp_rs(:, :, 3);
mask = (r_temp == 0) & (g_temp == 0) & (b_temp == 0);
r_temp(mask) = 255;
g_temp(mask) = 255;
b_temp(mask) = 255;
img_interp_rs(:, :, 1) = r_temp;
img_interp_rs(:, :, 2) = g_temp;
img_interp_rs(:, :, 3) = b_temp;

%Chromacity plot
sz = 11;
f1 = n_cols1/new_x_max;
f2 = n_rows1/new_y_max;

%Chart customization
chart_x_min_real = 0;
chart_x_max_real = new_x_max;
chart_y_min_real = 0;
chart_y_max_real = new_y_max;

chart_x_min = round(chart_x_min_real*f1);%unit conversion
chart_x_max = round(chart_x_max_real*f1);
chart_y_min = round(chart_y_min_real*f2);
chart_y_max = round(chart_y_max_real*f2);
%Tick marks
real_x = 0:.1:.8; 
real_y = 0:.1:.9;
x_temp_labels = round(real_x*f1);
y_temp_labels = round(real_y*f2);

%text
x_str = strsplit(num2str(real_x));
y_str = strsplit(num2str(real_y));

figure
imshow(img_interp_rs)
hold on
% plot(x_posi2*f1, y_posi2*f2, '.', 'Color', [1, 0, 1])
% plot(x_nega*f1, y_nega*f2, '.', 'Color', [0.9100 0.4100 0.1700])
plot(x_colour*f1, y_colour*f2, '-k', 'LineWidth', 2)
plot(xy_white(1)*f1, xy_white(2)*f2, '.', 'MarkerFaceColor', 'white', 'MarkerSize', sz)
plot(xy_R(1)*f1, xy_R(2)*f2, 'r.', 'MarkerSize', sz)
plot(xy_G(1)*f1, xy_G(2)*f2, 'g.', 'MarkerSize', sz)
plot(xy_B(1)*f1, xy_B(2)*f2, 'b.', 'MarkerSize', sz)
hold off
daspect([1 1 1])
ax = gca;
set(ax,'YDir','normal')
grid(ax, "on")
ax.GridColor = 'k';
ax.GridAlpha = .5; 
axis on
xlim([chart_x_min, chart_x_max])
ylim([chart_y_min, chart_y_max])
xticks(x_temp_labels);
yticks(y_temp_labels)
xticklabels(x_str)
yticklabels(y_str)

%Formatting
% set(zoom(gcf),'Motion','horizontal','Enable','on');
xlabel('x colour')
ylabel('y colour')
title('CIE chromaticity diagram')

%% Plot 1: Sensitivity function (comparing estimated varlues)

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

%% Plot 2: Spectral transmission interferences

chart_y2_min_real = 400;
chart_y2_max_real = 700;
chart_y2_range = chart_y2_max_real - chart_y2_min_real;

%coloured stripe
idx_temp = (lambda_var >= 400) & (lambda_var <= 700);
n_colors = sum(idx_temp);
lambda_var2 = lambda_var(idx_temp);
RGB_output = spectrumRGB(lambda_var2); %sRGB format
HSV_output = rgb2hsv(RGB_output);

%images
H_image = repmat(HSV_output(1, :, 1)', [1, 1001]);
sim_S_palette = repmat(0.5*ones([n_colors, 1]), [1, 1001]);
st_rescaled = rescale(spectral_transmission(idx_temp, :), 0, 1);
hsv_fig = cat(3, H_image, sim_S_palette, st_rescaled);
rgb_fig = hsv2rgb(hsv_fig);

st_chart_aspectR = (max_retardation-min_retardation)/(max_lambda-min_lambda);

st_rgb = rgb_fig;
% st_rgb = uint8(rescale(spectral_transmission, 0, 255));

numrows = size(st_rgb, 1);
numcols = round(st_chart_aspectR*numrows);
st_rs = imresize(st_rgb, [numrows, numcols]);

new_x_max = numcols;
new_y_max = n_colors;

%Tick marks
real_x = 0:500:max_retardation; 
real_y = chart_y2_min_real:100:chart_y2_max_real;
x_temp_labels = round(real_x*(new_x_max/max_retardation));
y_temp_labels = round(((real_y-400)*(new_y_max/chart_y2_range)));

%text
x_str = strsplit(num2str(real_x));
y_str = strsplit(num2str(real_y));

%Plot
hFig = figure;
imshow(st_rs)

%axis settings
ax = gca;
% set(ax,'YDir', 'normal')
grid(ax, "on")
ax.GridColor = 'k';
ax.GridAlpha = 1; 
axis on
xticks(x_temp_labels);
yticks(y_temp_labels)
xticklabels(x_str)
yticklabels(y_str)

%Formatting
% set(zoom(gcf),'Motion','horizontal','Enable','on');
xlabel('Retardation \Gamma (nm)')
ylabel('Light spectrum \lambda (nm)')
title('Spectral transmission interferences')

%Re-fitting within figure
set(gca, 'Position', [0.15,0.1, .7, 0.8])
set(gcf, 'OuterPosition', [100, 100, 1200, 400]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%% Plots 3: Michel-Levy

%x
new_x_min = 0;
new_x_max = size(rgb, 2); %range_fraction + 1
x_temp = new_x_min:new_x_max;
%y
new_y_max = round(range_fraction/chartAspectRatio); %true size
rgb_fig = repmat(rgb, [new_y_max, 1]); %image
[numrows, numcols, ~] = size(rgb_fig);
%stretching
extracted_fraction = (chart_y_max_real-chart_y_min_real)/(max_thickness-0); 
numrows = round(numrows/extracted_fraction);
rgb_fig = imresize(rgb_fig, [numrows, numcols]);

%Chart customization
chart_x_min = round(chart_x_min_real*(new_x_max/max_retardation));%unit conversion
chart_x_max = round(chart_x_max_real*(new_x_max/max_retardation));
chart_y_min = round(chart_y_min_real*(new_y_max/max_thickness));
chart_y_max = round(chart_y_max_real*(new_y_max/max_thickness));
%Tick marks
real_x = 0:500:max_retardation; 
real_y = 0:5:max_thickness;
x_temp_labels = round(real_x*(new_x_max/max_retardation));
y_temp_labels = round((real_y*(new_y_max/max_thickness))/extracted_fraction);

lines_bf = [0.0001, 0.005, 0.01, 0.015, 0.02, ...
    0.025, 0.03, 0.035, 0.04, 0.045, 0.05, ...
    0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2];
n_lines = length(lines_bf);
%prevent Inf (0 <--> 0.0001)

%text
x_str = strsplit(num2str(real_x));
y_str = strsplit(num2str(real_y));
lines_str = strsplit(num2str(lines_bf));

%extracted portion
stretched_y_min = chart_y_min/extracted_fraction;
stretched_y_max = chart_y_max/extracted_fraction;

hFig = figure;

imshow(rgb_fig, 'YData', [stretched_y_min, stretched_y_max])
hold on
text_x = zeros(1, n_lines);
text_y = zeros(1, n_lines);
for k = 1:n_lines
    f_temp = (x_temp/(lines_bf(k)*(chartAspectRatio/0.05)))/extracted_fraction;    
    plot(x_temp, f_temp, 'k-');
    
    %text locations
    if lines_bf(k) <= 0.05 %horizontal axis

        [x_data, y_data] = prepareCurveData(f_temp, x_temp);
        [fitresult, ~] = fit(x_data, y_data, 'linearinterp', 'Normalize', 'on');
        text_x(k) = fitresult(stretched_y_max);
        text_y(k) = stretched_y_max;

    else %vertical axis

        [x_data, y_data] = prepareCurveData(x_temp, f_temp);
        [fitresult, ~] = fit(x_data, y_data, 'linearinterp', 'Normalize', 'on');
        text_y(k) = fitresult(new_x_max);
        text_x(k) = new_x_max;
    end
    
end
t = text(text_x' + 5, text_y' + 5, lines_str');
hold off

%axis settings
ax = gca;
set(ax,'YDir', 'normal')
grid(ax, "on")
ax.GridColor = 'k';
ax.GridAlpha = 1; 
axis on
xlim([chart_x_min, chart_x_max])
ylim([stretched_y_min, stretched_y_max])
xticks(x_temp_labels);
yticks(y_temp_labels)
xticklabels(x_str)
yticklabels(y_str)

%Formatting
% set(zoom(gcf),'Motion','horizontal','Enable','on');
xlabel('Path difference (nm)')
ylabel('Sample thickness (microns)')
title('Calculated Michel-Levy Colour Chart')
subtitle('Birrefringence')

%Re-fitting within figure
set(gca, 'Position', [0.15,0.1, .7, 0.8])
set(gcf, 'OuterPosition', [100, 100, 1200, 600]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%% Plots 4: Raith–Sørensen chart

[fig_rows, fig_cols, n_channels] = size(rgb_fig);
n_pixels = fig_rows*fig_cols;
[row, col] = ind2sub([fig_rows, fig_cols], 1:n_pixels);
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

%sampling grid density
figure, 
scatter(trans_coord(:, 1), trans_coord(:, 2), 'r.')
xlim([0, val_top_x])
ylim([0, val_top_y])
grid on
axis equal

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

%Plot
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
xlim([chart_x_min, chart_x_max])
ylim([stretched_y_min, stretched_y_max])
xticks(x_temp_labels);
yticks(y_temp_labels)
xticklabels(x_str)
yticklabels(y_str)
% set(zoom(gcf), 'Motion', 'horizontal', 'Enable','on');
xlabel('Birrefringence')
ylabel('Sample thickness (microns)')
title('Raith–Sørensen chart')

%Re-fitting within figure
set(gca, 'Position', [0.15,0.1, .7, 0.8])
set(gcf, 'OuterPosition',[100, 100, 1200, 600]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
