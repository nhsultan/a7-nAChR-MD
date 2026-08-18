%% plot_M2_residues.m
% Description: Generates a tiled scatter plot of M2 lining residues, overlaid 
% with the HOLE pore radius profiles from both restrained and unrestrained runs.

% 1. Global Figure Settings
Size = 22;
figure(1)
t = tiledlayout(2, 2);
t.TileSpacing = 'none';
t.Padding = 'compact';
ylabel(t, "Channel Axis (\AA)", "FontSize", Size, "Interpreter", "latex", "FontName", "Times New Roman");
xlabel(t, "Radial Axis (\AA)", "FontSize", Size, "Interpreter", "latex", "FontName", "Times New Roman");

% 2. System and Data Definitions
resname = ["GLU", "SER", "THR", "LEU", "VAL", "LEU", "GLU"];
resid = [237, 240, 244, 247, 251, 254, 258];
folder = ["7KOX", "9LH5", "7EKT", "8V80"];
dim = ["X", "Y", "Z"];
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];

% Alignment offsets for E237 reference
SIM_E237 = [-22, -4, -20, -20];
PDB_E237 = [-22, -4, -4, -22];

axs = gobjects(1, length(folder));

% 3. Plotting Loop
for i = 1:length(folder)
    axs(i) = nexttile;
    hold on; box on;
    
    ax = gca;
    ax.FontSize = Size;
    ax.FontName = "Times New Roman";

    for ii = 1:length(resname)
        % Load X, Y, Z coordinate data
        for iii = 1:3
            FileName = append("../results/", folder(i), '/', folder(i), '.', resname(ii), '.', string(resid(ii)), '.', dim(iii), '.dat');
            T = table2array(readtable(FileName));

            % Convert coordinates to Angstroms
            if dim(iii) == "X"
                X = T .* 10;
            elseif dim(iii) == "Y"
                Y = T .* 10;
            else
                Z = T .* 10;
            end
        end

        % Calculate Mean Z and Radial Distances
        AvgZ = mean(Z, 2);
        R = zeros(size(T,1), size(T,2));

        for iii = 1:size(T,2)
            R(:,iii) = sqrt(X(:,iii).^2 + Y(:,iii).^2);
        end
        AvgR = min(R, [], 2);
        
        % 4. System-Specific Alignments & Scatter Plotting
        if folder(i) == "8V80"
            scatter(AvgR, -AvgZ, Marker=".");
            ax.YTick = [];
        elseif folder(i) == "7EKT"
            scatter(AvgR, AvgZ, Marker=".");
        elseif folder(i) == "9LH5"
            scatter(AvgR, AvgZ - 47.75, Marker=".");
            ax.XTick = [];
            ax.YTick = [];
        elseif folder(i) == "7KOX"
            scatter(AvgR, AvgZ, Marker=".");
            ax.XTick = [];
        end
    end
    
    % Set limits for all tiles
    xlim([0 13]);
    ylim([-5 40]);
    
    % 5. Overlay HOLE Pore Radius Profiles
    % Plot Unrestrained HOLE Profile (Solid line)
    Radius = readmatrix(append("../results/Unrestrained/Radius", folder(i), ".txt"));
    Radius(:,2) = Radius(:,2) - SIM_E237(i);

    if folder(i) == "8V80"
        plot(flip(Radius(:,1)), Radius(:,2), 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
    else
        plot(Radius(:,1), Radius(:,2), 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
    end
    
    % Plot Restrained HOLE Profile (Dashed line)
    RadiusPDB = readmatrix(append("../results/Restrained/Radius", folder(i), ".txt"));
    RadiusPDB(:,2) = RadiusPDB(:,2) - PDB_E237(i);
    
    plot(RadiusPDB(:,1), RadiusPDB(:,2), 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 1.5);
    yline(0, 'Color', 'k', 'LineStyle', '--');
end