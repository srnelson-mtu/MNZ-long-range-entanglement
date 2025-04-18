%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization and Dashboard %
clc
clear
close all;  

% Define density of states
densityofstates = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];

% Set Maximum Simulation Time
Tmax = 90; % Normalized Over G11
G11 = 1;   % Normalized G11

% Define constraints for analyzing Omegas
Ostep = ( 0.005) * G11; % Iteration accuracy ( _ ) -> Smaller is better
Olimit = 1.5 * G11;

% Set Time Evolution and Pump Intensity %
O1_values = -Olimit * G11:Ostep:Olimit * G11; % pump intensity applied on QDa
O2_values = -Olimit * G11:Ostep:Olimit * G11; % pump intensity applied on QDb
conc_results = zeros(length(O1_values), length(O2_values)); % Matrix to Store Concurrence Results

% Time evolution
timein = [0, Tmax];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Form Sigma_y Pauli Spin Matrices %
sy = zeros(4);
for t = 1:4
    for tt = 1:4
        if (t == 1 && tt == 2) || (t == 2 && tt == 1)
            sy(t, tt) = 1;
        elseif (t == 3 && tt == 4) || (t == 4 && tt == 3)
            sy(t, tt) = -1;
        else
            sy(t, tt) = 0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve Master Equations and Calculate Concurrence %
for i = 1:length(O1_values)
    for j = 1:length(O2_values)
        % Set current O1 and O2 values
        O1 = O1_values(i);
        O2 = O2_values(j);
        
        % Solve the master equation
        [T, Y] = ode15s(@(t, rho) master_equation_test(t, rho, O1, O2), [0, Tmax], densityofstates);
        
        % Calculate concurrence
        conc = zeros(size(T));
        for time = 1:length(T)
            Rho = reshape(Y(time, :), 4, 4); 
            Rho_tilda = sy * conj(Rho) * sy;      % rho_tilda
            EV = abs(sqrt(eig(Rho * Rho_tilda))); % sqrt of the eign values of rho*rho_tilda
            EV = sort(EV);                        % sorting the sqrt of the eigen values of rho*rho_tilda
            conc(time) = max([0, EV(4) - EV(3) - EV(2) - EV(1)]); % C(t)
        end
        
        % Save the last concurrence value
        conc_results(i, j) = conc(end);
        
         % Extract the last state
        Rho_end = reshape(Y(end, :), 4, 4);
              
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the maximum concurrence value
[max_concurrence, max_index] = max(conc_results, [], 'all');

% Find the row and column indices corresponding to the maximum value
[max_row, max_col] = ind2sub(size(conc_results), max_index);

% Retrieve the Omega values corresponding to the maximum concurrence
O1_max_concurrence = O1_values(max_row);
O2_max_concurrence = O2_values(max_col);

% Display the maximum concurrence value
disp(['Maximum Concurrence Value: ', num2str(max_concurrence)]);

% Display the Omega values for maximum concurrence
disp(['Omega values (O1 and O2) for maximum concurrence:']);
disp(['O1 = ', num2str(O1_max_concurrence)]);
disp(['O2 = ', num2str(O2_max_concurrence)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Concurrence Results %
figure;

% Plot the concurrence results as a heatmap
imagesc(O1_values, O2_values, conc_results);
colormap(parula); % Plot the concurrence as a custom colormap

% Add color bar with custom font size
h = colorbar;
h.FontSize = 18; % Set color bar font size
title(h, 'C_{ss}'); % Add color bar title

% Set color bar limits
caxis([0, 0.4]); % Set the color bar limits (from 0 to 0.4)
h.Ticks = 0:0.1:0.4; % Set color bar tick increments (0, 0.1, 0.2, 0.3, 0.4)

% Adjust font properties for axes and title
xlabel('\Omega_1/\gamma', 'FontName', 'Arial', 'FontSize', 16);
ylabel('\Omega_2/\gamma', 'FontName', 'Arial', 'FontSize', 16);
title('r_{12}/\lambda_0=0.5', 'FontName', 'Arial', 'FontSize', 20); % Change Title HERE 

% Adjust font properties for tick labels
ax = gca;
ax.FontName = 'Arial'; % Set axes font to Arial
ax.FontSize = 16; % Set axes font size
ax.XTick = -1.5:0.5:1.5; % Set X-axis tick increments (-1.5, 0, 1.5)
ax.YTick = -1.5:0.5:1.5; % Set Y-axis tick increments (-1.5, 0, 1.5)

% Adjust font properties for color bar tick labels
h.FontName = 'Arial'; % Set color bar font to Arial

% Add figure identifier "(a)" to the upper left corner with adjusted position
text(ax.XLim(1), ax.YLim(2) + 0.05*(ax.YLim(2)-ax.YLim(1)), '(a)', ... 
    'FontSize', 16, 'FontWeight', 'bold','FontName', 'Arial');

% Set aspect ratio to ensure proper display
axis xy;

% Ensure figure size is appropriate
set(gcf, 'Position', [100, 100, 800, 600]); % Set figure size

% Save figure if needed
saveas(gcf, 'concurrence_figure.png'); % Save figure as PNG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Master Equation Function &
function drho = master_equation_test(t, rho, O1, O2)
    % Constants
    eps0 = 8.854e-12;
    deb = 3.33564e-30;
    hbar = 1.054e-34;
    d1 = 60;
    d2 = 60;
    
    % Dipole moments in Debye
    dip1 = d1 * deb;
    dip2 = d2 * deb;
    
    % Define coupling constants % Change Variables HERE -> 
    G11 = 2 * dip1 * dip1 / (eps0 * hbar);
    G22 = 2 * dip1 * dip1 / (G11 * eps0 * hbar);
    G12 = 2 * dip1 * dip2 * ( 0.95 ) / (G11 * eps0 * hbar); % Change ( _ ) Number, Previous Program: 0.8
    % For r12/lambda0 = 0.5  -> New: 0.95   Old: 0.85
    % For r12/lambda0 = 1.0  -> New: 0.9    Old: 0.50
    % For r12/lambda0 = 1.5  -> New: 0.85   Old: 0.10
    % For r12/lambda0 = 5.0  -> New: 0.75   Old:
    % For r12/lambda0 = 10.0 -> New: 0.5    Old:
    % For r12/lambda0 = 15.0 -> New: 0.15   Old:
    G21 = G12;
    g12 = dip1 * dip2 * ( -0.70 ) / (G11 * eps0 * hbar); % Change ( _ ) Number, Previous Program: 0.5
    % For r12/lambda0 = 0.5  -> New: -0.70  Old: -0.15
    % For r12/lambda0 = 1.0  -> New: -0.71  Old: -0.20
    % For r12/lambda0 = 1.5  -> New: -0.72  Old: -0.25
    % For r12/lambda0 = 5.0  -> New: -0.71  Old:
    % For r12/lambda0 = 10.0 -> New: -0.65  Old:
    % For r12/lambda0 = 15.0 -> New: -0.5   Old:
    g21 = g12;
    
    % Redefine G11 to Renormalize the System
    G11 = 1;
    
    % Master equation
    drho = zeros(16, 1);
    drho(1) = 1i*O1*(rho(13)-rho(4)) + 1i*O2*(rho(9)-rho(3)) + G11*rho(16) + G22*rho(11) + G12*rho(12) + G21*rho(15);
    drho(2) = 1i*O1*(rho(14)-rho(3)) + 1i*O2*(rho(10)-rho(4)) - 0.5*rho(2)*(G11+G22);
    drho(3) = 1i*O1*(rho(15)-rho(2)) + 1i*O2*(rho(11)-rho(1)) + G11*rho(14) - 0.5*G22*rho(3) - 0.5*G12*rho(4) + G12*rho(10) - 1i*g12*rho(4);
    drho(4) = 1i*O1*(rho(16)-rho(1)) + 1i*O2*(rho(12)-rho(2)) - 0.5*G11*rho(4) + G22*rho(10) - 0.5*G21*rho(3) + G21*rho(14) - 1i*g21*rho(3);
    
    drho(5) = -1i*O1*(rho(8)-rho(9)) - 1i*O2*(rho(7)-rho(13)) - 0.5*rho(5)*(G11+G22);
    drho(6) = -1i*O1*(rho(7)-rho(10)) - 1i*O2*(rho(8)-rho(14)) - rho(6)*(G11+G22);
    drho(7) = 1i*O1*(rho(11)-rho(6)) - 1i*O2*(rho(5)-rho(15)) - rho(7)*(0.5*G11+G22) - 0.5*G12*rho(8) - 1i*g12*rho(8);
    drho(8) = -1i*O1*(rho(5)-rho(12)) + 1i*O2*(rho(16)-rho(6)) - rho(8)*(G11+0.5*G22) - 0.5*G21*rho(7) - 1i*g21*rho(7);
    
    drho(9) = 1i*O1*(rho(5)-rho(12)) + 1i*O2*(rho(1)-rho(11)) + G11*rho(8) - 0.5*G22*rho(9) + G21*rho(7) - 0.5*G21*rho(13) + 1i*g21*rho(13);
    drho(10) = 1i*O1*(rho(6)-rho(11)) - 1i*O2*(rho(12)-rho(2)) - rho(10)*(0.5*G11+G22) - 0.5*G21*rho(14) + 1i*g21*rho(14);
    drho(11) = 1i*O1*(rho(7)-rho(10))-1i*O2*(rho(9)-rho(3)) + G11*rho(6) - G22*rho(11) - 0.5*G12*rho(12) - 0.5*G21*rho(15) + 1i*(g21*rho(15)-g12*rho(12));
    drho(12) = 1i*O1*(rho(8)-rho(9)) + 1i*O2*(rho(4)-rho(10)) - 0.5*rho(12)*(G11+G22) + G21*rho(6) - 0.5*G21*(rho(11)+rho(16)) + 1i*g21*(rho(16)-rho(11));
    
    drho(13) = 1i*O1*(rho(1)-rho(16)) + 1i*O2*(rho(5)-rho(15)) - 0.5*G11*rho(13) + G22*rho(7) - 0.5*G12*rho(9) + G12*rho(8) + 1i*g12*rho(9);
    drho(14) = -1i*O1*(rho(15)-rho(2)) + 1i*O2*(rho(6)-rho(16)) - rho(14)*(G11+0.5*G22) - 0.5*G12*rho(10) + 1i*g12*rho(10);
    drho(15) = 1i*O1*(rho(3)-rho(14)) + 1i*O2*(rho(7)-rho(13)) - 0.5*rho(15)*(G11+G22) + 0.5*G12*(2*rho(6)-rho(11)-rho(16)) + 1i*g12*(rho(11)-rho(16));
    drho(16) = -1i*O1*(rho(13)-rho(4)) + 1i*O2*(rho(8)-rho(14)) -G11*rho(16) + G22*rho(6) - 0.5*G12*rho(12) - 0.5*G21*rho(15) + 1i*(g12*rho(12)-g21*rho(15));
 
end
