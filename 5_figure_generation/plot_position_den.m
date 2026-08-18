%% plot_position_den.m
% Description: Generates the 1x5 tiled layout Position Density plot for 
% specific residues along the channel axis of the a7-nAChR.

% 1. Global Figure Settings
Size = 22;
linewidth = 1.5;
fontsize = 18;

% Create figure and tiled layout
figure('Position', [744 243 560 835])
t = tiledlayout(1, 5);
xlabel(t, "Position Density, P$_{r}$ (-)", "FontSize", Size, "Interpreter", "latex", "FontName", "Times New Roman");
t.TileSpacing = 'none';
t.Padding = 'compact';

% 2. System and Residue Definitions
ResID = ["41", "43", "81", "88", "96", "97", "98", "100", "104", "237", "258"];
ResNameID = ["D:41", "D:43", "D:81", "D:88", "D:96", "E:97", "R:98", "D:100", "H:104", "E:237", "E:258"];
ResIDPosition = zeros(length(ResID), 2);
PDB = ["7KOX", "7EKT", "8V80", "8V82", "9LH5"];
color = ['b', 'r', 'g', 'c', 'm'];

% Z-axis Alignment Offsets
Diff = [0, 0, 0, 0, -80];
E237 = [49.75, 47.75, 61.75, 52.75, 140.25];

% Array to store axes for linking later
axs = gobjects(1, length(PDB));

% 3. Plotting Loop
for ii = 1:length(PDB)
    axs(ii) = nexttile;
    hold on; box on;

    for i = 1:length(ResID)
        % Construct dynamic filename based on pipeline output
        filename = append("../results/", PDB(ii), "/resDensity.", PDB(ii), '.', ResID(i), ".sidechain.dat");

        % Import and normalize the density data
        resDensity = readmatrix(filename);
        resDensity(:,2) = resDensity(:,2) ./ sum(resDensity(:,2));

        % Find peak density for labeling
        [yMax, maxIdx] = max(resDensity(:,2));
        xMax = resDensity(maxIdx, 1);

        % Apply alignment offsets
        ResIDPosition(i,1) = xMax - Diff(ii) - E237(ii);
        ResIDPosition(i,2) = yMax;
        resDensity(:,1) = resDensity(:,1) - Diff(ii) - E237(ii);
        
        % Plot density curve
        plot(resDensity(:,2), resDensity(:,1), 'Color', color(ii), 'LineStyle', '-', 'LineWidth', linewidth);
    end
    
    % Add residue text labels at peak positions
    text(ResIDPosition(:,2), ResIDPosition(:,1), ResNameID, "Interpreter", "latex", "FontName", "Times New Roman", "FontSize", fontsize, "Color", 'k');

    % 4. Draw Shaded Region (Between E237 and E258)
    AreaShadeX = [0 0 0.3 0.3];
    % Indices 10 and 11 correspond to E237 and E258
    AreaShadeY = [ResIDPosition(10,1) ResIDPosition(11,1) ResIDPosition(11,1) ResIDPosition(10,1)];
    patch(AreaShadeX, AreaShadeY, color(ii), 'FaceAlpha', 0.15, 'LineStyle', ':');

    % 5. Format Axes
    ylim([-15 95]);
    xlim([0 0.30]);
    grid off;
    xline(0, 'Color', color(ii), 'LineStyle', '-', 'LineWidth', linewidth);
    
    ax = gca;
    ax.FontSize = Size;
    ax.FontName = "Times New Roman";
    
    if ii == 5
        set(gca, 'YDir', 'reverse');
    end
    
    if ii == 1
        ylabel("Channel Axis (\AA)", "FontSize", Size, "Interpreter", "latex", "FontName", "Times New Roman");
    else
        ax.YTick = []; % Hide Y-ticks for inner panels
    end
end

% 6. Link all Y-axes for unified zooming/panning
linkaxes(axs, 'y');