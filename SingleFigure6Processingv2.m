% List of .fig files to be combined
figFiles = {'r12lambda_0.5.fig', 'r12lambda_1.0.fig', ...
            'r12lambda_1.5.fig', 'r12lambda_5.0.fig', ...
            'r12lambda_10.0.fig', 'r12lambda_15.0.fig'};

% Number of rows and columns for the subplot layout
numRows = 2;
numCols = 3;

% Create a new figure for the combined output
combinedFig = figure;
set(combinedFig, 'Position', [100, 100, 1600, 800]);  % Slightly smaller figure size

% Initialize global color limits
cMin = inf;
cMax = -inf;

% Editable text properties
fontName = 'Arial';      
fontSize = 22;           
titleFontSize = 24;      
labelFontSize = 22;      
annotationFontSize = 22; 
colorBarFontSize = 22;   
tickLabelFontSize = 18;  
colorBarTickFontSize = 18;

% Extract data from each .fig and store it in a structure
axesData = struct();
for i = 1:length(figFiles)
    try
        % Open the figure invisibly
        fig = openfig(figFiles{i}, 'invisible');
        
        % Find axes and image objects in the figure
        ax = findobj(fig, 'Type', 'axes');
        imgObj = findobj(ax, 'Type', 'image');
        
        % Extract the heatmap data
        XData = get(imgObj, 'XData');
        YData = get(imgObj, 'YData');
        CData = get(imgObj, 'CData');
        
        % Update global color limits
        cMin = min(cMin, min(CData(:)));
        cMax = max(cMax, max(CData(:)));
        
        % Store the axes handle and data
        axesData(i).Title = get(get(ax, 'Title'), 'String');
        axesData(i).TitleInterpreter = get(get(ax, 'Title'), 'Interpreter');
        axesData(i).XLabel = get(get(ax, 'XLabel'), 'String');
        axesData(i).XLabelInterpreter = get(get(ax, 'XLabel'), 'Interpreter');
        axesData(i).YLabel = get(get(ax, 'YLabel'), 'String');
        axesData(i).YLabelInterpreter = get(get(ax, 'YLabel'), 'Interpreter');
        axesData(i).XData = XData;
        axesData(i).YData = YData;
        axesData(i).CData = CData;
        axesData(i).Annotations = findobj(ax, 'Type', 'text');
        
        % Close the original figure
        close(fig);
    catch ME
        warning('Error processing %s: %s', figFiles{i}, ME.message);
    end
end

% Create a tiled layout for better control of spacing
t = tiledlayout(numRows, numCols, 'TileSpacing', 'compact', 'Padding', 'compact');

% Adjust the outer position to add more spacing
t.OuterPosition = [0.02, 0.05, 0.85, 0.88];  % Slightly smaller usable space for subplots

% Create subplots for each figure
for i = 1:length(axesData)
    try
        % Skip if no data available
        if isempty(axesData(i).CData)
            warning('Skipping empty data for figure %d.', i);
            continue;
        end
        
        % Create a subplot in the tiled layout
        nexttile;
        
        % Plot heatmap
        imagesc(axesData(i).XData, axesData(i).YData, axesData(i).CData);
        axis xy tight;
        
        % Set custom color bar limits
        caxis([0, 0.4]);  % Set color bar limits to range from 0 to 0.4
        
        % Add title, labels, and formatting
        title(axesData(i).Title, 'Interpreter', axesData(i).TitleInterpreter, 'FontSize', titleFontSize, 'FontName', fontName);
        xlabel(axesData(i).XLabel, 'Interpreter', axesData(i).XLabelInterpreter, 'FontSize', labelFontSize, 'FontName', fontName);
        ylabel(axesData(i).YLabel, 'Interpreter', axesData(i).YLabelInterpreter, 'FontSize', labelFontSize, 'FontName', fontName);
        
        % Adjust tick label font size and style
        ax = gca;  % Get current axes handle
        set(ax, 'FontSize', tickLabelFontSize, 'FontName', fontName);
        
        % Set ticks to -1, 0, and 1 only
        ax.XTick = [-1, 0, 1];
        ax.YTick = [-1, 0, 1];
        
        % Ensure tick labels are horizontal and simple
        ax.XTickLabelRotation = 0;  % No rotation
        ax.YTickLabelRotation = 0;  % No rotation

        % Add notation (a), (b), etc. with bold font
        annotationText = sprintf('(%c)', 'a' + i - 1); % Generate (a), (b), etc.
        text(-1.5, 1.7, annotationText, 'FontSize', annotationFontSize, 'FontWeight', 'bold', 'FontName', fontName, 'Interpreter', 'latex');
        
        % Re-add annotations with preserved positions
        for annotation = axesData(i).Annotations'
            annotationString = get(annotation, 'String');
            annotationPos = get(annotation, 'Position');
            text(annotationPos(1), annotationPos(2), annotationString, ...
                'FontSize', get(annotation, 'FontSize'), ...
                'FontWeight', get(annotation, 'FontWeight'), ...
                'FontName', get(annotation, 'FontName'), ...
                'Interpreter', get(annotation, 'Interpreter'));
        end
    catch ME
        warning('Error creating subplot for figure %d: %s', i, ME.message);
    end
end

% Create the color bar
cb = colorbar('eastoutside');
cb.Position = [0.88, 0.1, 0.02, 0.8];  % Adjust position for proper spacing

% Set color bar ticks and limits
cb.Ticks = 0:0.1:0.4;
set(cb, 'FontSize', colorBarTickFontSize, 'FontName', fontName);

% Add the color bar title
axes('Position', [0.88, 0.92, 0.02, 0.05], 'Visible', 'off');  % Invisible axis for the title
text(0.5, 0.5, '$C_{ss}$', ...
    'Interpreter', 'latex', 'FontSize', colorBarFontSize, 'FontName', fontName, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
