%% plot_channel_radius_comparison.m
% Description: 
% Plots the average channel radius data for five PDB models across three 
% different conditions (Restrained, Unrestrained, Modulator) in a 1x3 tiled layout.
%
% DATA INTERPRETATION:
% - Column 1: Channel Radius (Angstroms)
% - Column 2: Axis Length / Position along the channel (Angstroms)

clc; close all; clear all;

%% 1. Configuration
pdb_ids = {'7EKT', '7KOX', '8V80', '8V82', '9LH5'}; 
system_folders = {'Restrained', 'Unrestrained', 'Modulator'}; 

% Line properties: Blue, Red, Green, Magenta, Cyan
line_colors = {'blue', 'red', 'green', 'magenta', 'cyan'}; 
line_styles = {'-', '-', '-', '-', '-'}; 

x_label_text = 'Channel Radius (\AA)'; 
y_label_text = 'Channel Axis (\AA)'; 

axis_handles = gobjects(1, length(system_folders)); 

%% 2. Initialize Figure and Tiled Layout
f = figure('Position', [50, 100, 1400, 500]); 
t = tiledlayout(1, 3, 'TileSpacing', 'none', 'Padding', 'compact'); 

%% 3. Loop through Systems and Plot Data
for sys_idx = 1:length(system_folders)
    current_folder = system_folders{sys_idx}; 

    ax = nexttile; 
    axis_handles(sys_idx) = ax; 
    hold on; 
    
    % Add translucent color patch to highlight the pore region
    patch_color = [0.6 0.6 0.6]; 
    patch([0 10 10 0], [0 0 32 32], patch_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5); 

    line_handles_for_legend = gobjects(0); 
    
    for pdb_idx = 1:length(pdb_ids)
        pdb = pdb_ids{pdb_idx}; 
        % Assuming data is inside the structured results folder
        filepath = fullfile('../results', current_folder, ['Radius.', pdb, '.txt']); 

        if exist(filepath, 'file')
            data = readmatrix(filepath); 
            cleaned_data = data(~isnan(data(:, 1)), :); 
            
            if size(cleaned_data, 1) > 0 
                radius = cleaned_data(:, 1); 
                x_axis = cleaned_data(:, 2); 
                
                % =========================================================
                % SPECIFIC DATA SHIFTS/FLIPS FOR STRUCTURAL ALIGNMENT
                % =========================================================
                
                % Flip axis length for Unrestrained 8V80
                if strcmp(pdb, '8V80') && strcmp(current_folder, 'Unrestrained') 
                    x_axis = flipud(x_axis); 
                end 

                % Alignment shifts for 9LH5
                if strcmp(pdb, '9LH5') 
                     if strcmp(current_folder, 'Unrestrained') 
                        x_axis = x_axis - 17; 
                     elseif strcmp(current_folder, 'Restrained') 
                        x_axis = x_axis - 16; 
                     elseif strcmp(current_folder, 'Modulator') 
                        x_axis = x_axis + 2; 
                     end 
                end 

                % Alignment shifts for 7EKT
                if strcmp(pdb, '7EKT') && (strcmp(current_folder, 'Restrained') || strcmp(current_folder, 'Modulator')) 
                    x_axis = x_axis - 17; 
                end 

                % Alignment shifts for 8V80 & 8V82
                if (strcmp(pdb, '8V80') || strcmp(pdb, '8V82')) && strcmp(current_folder, 'Restrained') 
                    x_axis = x_axis + 2; 
                end 
                
                % Global offset to normalize to E-237
                x_axis = x_axis + 20; 

                % =========================================================

                % Plot the profile line
                h_line = plot(radius, x_axis, ...
                    'DisplayName', pdb, ...
                    'Color', line_colors{pdb_idx}, ...
                    'LineStyle', line_styles{pdb_idx}, ...
                    'LineWidth', 1.5); 
                
                line_handles_for_legend(end+1) = h_line; 
            end 
        end 
    end 
    
    % --- Subplot Customization ---
    xlabel(x_label_text, 'FontSize', 14, 'Interpreter', 'latex'); 
    grid off; box on; 
    ylim([-10 80]); xlim([0 10]); 

    % Panel Labels
    if sys_idx == 1 
        ylabel(y_label_text, 'FontSize', 14, 'Interpreter', 'latex'); 
        if ~isempty(line_handles_for_legend) 
            legend(line_handles_for_legend, pdb_ids, 'Location', 'best', 'Box', 'off');  
        end 
        text(0.5, 75, '(a) Restrained', 'FontSize', 14, 'Interpreter', 'latex'); 
            
    elseif sys_idx == 2 
        yticklabels([]); yticks([]); 
        text(0.5, 75, '(b) Unrestrained', 'FontSize', 14, 'Interpreter', 'latex'); 
        
    elseif sys_idx == 3 
        yticklabels([]); yticks([]); 
        text(0.5, 75, '(c) With Modulators', 'FontSize', 14, 'Interpreter', 'latex'); 
    end 
    hold off; 
end 

linkaxes(axis_handles, 'y'); 