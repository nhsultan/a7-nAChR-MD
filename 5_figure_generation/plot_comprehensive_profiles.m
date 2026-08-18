%% plot_comprehensive_profiles.m
% Description: Generates a 3x4 grid plotting Electrostatics, Ion Density, 
% Water Density, and Hydration Numbers.
% Rows: 1=Restrained, 2=Unrestrained, 3=Modulators.

clc; clear all; close all;

%% 1. Configuration & Shading Preprocessing
ResID = ["237", "258"]; % Residues used to define the shaded pore region
PDB = ["7KOX", "7EKT", "8V80", "8V82", "9LH5"]; 
color = ['b', 'r', 'g', 'c', 'm']; 
E237 = [49.75, 47.75, 61.75, 52.75, 110.25]; % Z-axis alignment offsets

Size = 22; 
fontsize = 18; 

Dirs = ["Restrained", "Unrestrained", "Modulator"]; 
AreaShades = cell(3, 5); % Cell array to store shading Y-coordinates for (Condition, PDB)

% Precalculate the shaded pore regions based on sidechain density peaks
for d = 1:length(Dirs) 
    for i = 1:length(PDB) 
        currentResPos = zeros(length(ResID), 2); 
        for ii = 1:length(ResID) 
            filename = fullfile("../results", Dirs(d), PDB(i), append("resDensity.", PDB(i), ".", ResID(ii), ".sidechain.dat")); 
            if exist(filename, 'file')
                resDensity = readmatrix(filename); 
                resDensity(:,2) = resDensity(:,2) ./ sum(resDensity(:,2)); 
                [yMax, maxIdx] = max(resDensity(:,2)); 
                currentResPos(ii,1) = resDensity(maxIdx,1) - E237(i); 
            end
        end
        % Shade boundary between GLU 237 and GLU 258
        AreaShades{d, i} = [currentResPos(1,1) currentResPos(2,1) currentResPos(2,1) currentResPos(1,1)]; 
    end 
end 

%% 2. Figure Setup
figure('Position', [50 50 1600 1200]) 
tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'tight') 

% Y-axis shift mappings per PDB based on the original script
% Format: [7KOX, 7EKT, 8V80, 8V82, 9LH5]
shift_ES_R = [18, 0, 22, 19, 0]; 
shift_ID_R = [20, 4, 20, 30, 15]; 
shift_WD_R = [25, 25, 25, 28, 25]; 
shift_HY_R = [20, 0, 10, 20, 5]; 

shift_ES_U = [20, 15, 20, 17, 0]; 
shift_ID_U = [22, 18, 22, 22, 5]; 
shift_WD_U = [20, 20, 35, 20, 5]; 
shift_HY_U = [0, 0, 0, 15, 0]; 

shift_ES_M = [20, 0, 20, 20, 5]; 
shift_ID_M = [2, 20, 0, 20, 20]; 
shift_WD_M = [22, 2, 20, 22, 6]; 
shift_HY_M = [15, 0, 15, 15, 0]; 

%% 3. Plotting Loop (Iterate over Rows/Conditions)
for row = 1:3
    current_dir = Dirs(row);
    
    % --- Column 1: Electrostatics (PME) ---
    nexttile; hold on; box on; 
    for i = 1:length(PDB) 
        filename = fullfile("../results", current_dir, PDB(i), append("pme.", PDB(i), ".dat")); 
        if exist(filename, 'file')
            pme = readmatrix(filename); x = -80:79; 
            
            % Select correct shift array
            if row==1; shift = shift_ES_R(i); elseif row==2; shift = shift_ES_U(i); else; shift = shift_ES_M(i); end 
            
            % Plotting (using column 2 for Modulators, column 1 for others, with 8V80 flipped in rows 1/2)
            if row == 3 || PDB(i) == "8V82" 
                plot(pme(:,2), x + shift, 'Color', color(i), 'LineWidth', 1.5); 
            elseif PDB(i) == "8V80" 
                plot(flip(pme(:,1)), x + shift, 'Color', color(i), 'LineWidth', 1.5); 
            else 
                plot(pme(:,1), x + shift, 'Color', color(i), 'LineWidth', 1.5); 
            end 
            patch([-10 -10 22 22], AreaShades{row, i}, color(i), 'FaceAlpha', 0.15, 'LineStyle', ':'); 
        end
    end
    ylim([-15 75]); xlim([-10 16]); 
    if row == 3; xlim([-10 22]); xlabel("Electric Potential (k$_{B}$T/$e^-$)", "FontSize", Size, "Interpreter", "latex"); end 
    ylabel("Channel Axis (\AA)", "FontSize", Size, "Interpreter", "latex"); 
    set(gca, 'FontSize', fontsize, 'FontName', 'Times New Roman'); 

    % --- Column 2: Ion Density ---
    nexttile; hold on; box on; 
    for i = 1:length(PDB) 
        filename = fullfile("../results", current_dir, PDB(i), append(PDB(i), ".pot.txt")); 
        if exist(filename, 'file')
            IonDensity = readmatrix(filename); 
            if row==1; den_scale = 4; else; den_scale = 20; end 
            Density = (IonDensity(:,2) / den_scale) / 0.10; 
            
            if row==1; shift = shift_ID_R(i); elseif row==2; shift = shift_ID_U(i); else; shift = shift_ID_M(i); end 
            
            % Flip logic for specific PDBs
            if (row == 1 || row == 2) && PDB(i) == "8V80" 
                plot(flip(Density), IonDensity(:,1) + shift, 'Color', color(i), 'LineWidth', 1.5); 
            elseif row == 3 && PDB(i) == "9LH5" 
                plot(Density, flip(IonDensity(:,1)) + shift, 'Color', color(i), 'LineWidth', 1.5); 
            else 
                plot(Density, IonDensity(:,1) + shift, 'Color', color(i), 'LineWidth', 1.5); 
            end 
            patch([0 0 22 22], AreaShades{row, i}, color(i), 'FaceAlpha', 0.15, 'LineStyle', ':'); 
        end
    end
    ylim([-15 75]); xlim([0 12]); set(gca, 'YTick', [], 'FontSize', fontsize, 'FontName', 'Times New Roman'); 
    if row == 3; xlabel("Ion Density, N$_{K^{+}}$ (-)", "FontSize", Size, "Interpreter", "latex"); end 

    % --- Column 3: Water Density ---
    nexttile; hold on; box on; 
    for i = 1:length(PDB) 
        filename = fullfile("../results", current_dir, PDB(i), append(PDB(i), ".water.txt")); 
        if exist(filename, 'file')
            WaterDensity = readmatrix(filename); 
            if row==1; water = movmean(WaterDensity(:,2), 5) ./ 12.5; else; water = WaterDensity(:,2) ./ 12.5; end 
            if row==1; shift = shift_WD_R(i); elseif row==2; shift = shift_WD_U(i); else; shift = shift_WD_M(i); end 
            
            if row == 2 && PDB(i) == "8V80" 
                plot(flip(water), WaterDensity(:,1) + shift, 'Color', color(i), 'LineWidth', 1.5); 
            else 
                plot(water, WaterDensity(:,1) + shift, 'Color', color(i), 'LineWidth', 1.5); 
            end 
            patch([0 0 22 22], AreaShades{row, i}, color(i), 'FaceAlpha', 0.15, 'LineStyle', ':'); 
        end
    end
    ylim([-15 75]); xlim([0 1]); xticks([0.2 0.4 0.6 0.8]); set(gca, 'YTick', [], 'FontSize', fontsize, 'FontName', 'Times New Roman'); 
    if row == 3; xlabel("Water Density, N$_{H_{2}O}$ (-)", "FontSize", Size, "Interpreter", "latex"); end 

    % --- Column 4: Hydration ---
    nexttile; hold on; box on; 
    HydrationTypes = ["protein.hydration", "water.hydration"]; 
    if row == 1 || row == 2; HydrationTypes = ["proteinhydration", "waterhydration"]; end 
    
    for i = 1:length(PDB) 
        % Filter based on your original script logic for Unrestrained
        if row == 2 && ~(PDB(i) == "9LH5" || PDB(i) == "7KOX" || PDB(i) == "8V82"); continue; end 
        
        for ii = 1:length(HydrationTypes) 
            filename = fullfile("../results", current_dir, PDB(i), append(PDB(i), ".", HydrationTypes(ii), ".txt")); 
            if exist(filename, 'file')
                hydration = readmatrix(filename); 
                if row==1; shift = shift_HY_R(i); elseif row==2; shift = shift_HY_U(i); else; shift = shift_HY_M(i); end 
                
                if row == 3 && PDB(i) ~= "9LH5" 
                    plot(hydration(:,2), flip(hydration(:,1)) + shift, 'Color', color(i), 'LineWidth', 1.5); 
                else 
                    plot(hydration(:,2), hydration(:,1) + shift, 'Color', color(i), 'LineWidth', 1.5); 
                end 
            end
        end
        patch([0 0 7 7], AreaShades{row, i}, color(i), 'FaceAlpha', 0.15, 'LineStyle', ':'); 
    end
    ylim([-15 75]); xlim([0 7]); xticks([1 3 5 7]); set(gca, 'YTick', [], 'FontSize', fontsize, 'FontName', 'Times New Roman'); 
    if row == 3; xlabel("Hydration Number (-)", "FontSize", Size, "Interpreter", "latex"); end 
end 