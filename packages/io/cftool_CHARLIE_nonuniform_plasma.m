function [fitresult, gof] = createFits(z_ne, ne_z, z_Te, Te_z, r_ne, ne_r, r_Te, Te_r)
%CREATEFITS(Z_NE,NE_Z,Z_TE,TE_Z,R_NE,NE_R,R_TE,TE_R)
%  Create fits.
%
%  Data for 'ne(z)-polyfit' fit:
%      X Input : z_ne
%      Y Output: ne_z
%  Data for 'Te(z)-Lorentzian' fit:
%      X Input : z_Te
%      Y Output: Te_z
%  Data for 'ne(z)-Lorentzian' fit:
%      X Input : z_ne
%      Y Output: ne_z
%  Data for 'ne(r)-Bessel' fit:
%      X Input : r_ne
%      Y Output: ne_r
%  Data for 'ne(r)-polyfit' fit:
%      X Input : r_ne
%      Y Output: ne_r
%  Data for 'Te(r)-polyfit' fit:
%      X Input : r_Te
%      Y Output: Te_r
%  Output:
%      fitresult : a cell-array of fit objects representing the fits.
%      gof : structure array with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 19-Apr-2021 19:59:58 自动生成

%% Initialization.

% Initialize arrays to store fits and goodness-of-fit.
fitresult = cell( 6, 1 );
gof = struct( 'sse', cell( 6, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

%% Fit: 'ne(z)-polyfit'.
[xData, yData] = prepareCurveData( z_ne, ne_z );

% Set up fittype and options.
ft = fittype( 'poly4' );

% Fit model to data.
[fitresult{1}, gof(1)] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', 'ne(z)-polyfit' );
h = plot( fitresult{1}, xData, yData );
legend( h, 'ne_z vs. z_ne', 'ne(z)-polyfit', 'Location', 'NorthEast' );
% Label axes
xlabel z_ne
ylabel ne_z
grid on

%% Fit: 'Te(z)-Lorentzian'.
[xData, yData] = prepareCurveData( z_Te, Te_z );

% Set up fittype and options.
ft = fittype( 'k*Lorentzian([0,half_gamma],x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.652420835422264 0.531048838345102];

% Fit model to data.
[fitresult{2}, gof(2)] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'Te(z)-Lorentzian' );
h = plot( fitresult{2}, xData, yData );
legend( h, 'Te_z vs. z_Te', 'Te(z)-Lorentzian', 'Location', 'NorthEast' );
% Label axes
xlabel z_Te
ylabel Te_z
grid on

%% Fit: 'ne(z)-Lorentzian'.
[xData, yData] = prepareCurveData( z_ne, ne_z );

% Set up fittype and options.
ft = fittype( 'k*Lorentzian([0,half_gamma],x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.935975553324428 0.389281682796525];

% Fit model to data.
[fitresult{3}, gof(3)] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'ne(z)-Lorentzian' );
h = plot( fitresult{3}, xData, yData );
legend( h, 'ne_z vs. z_ne', 'ne(z)-Lorentzian', 'Location', 'NorthEast' );
% Label axes
xlabel z_ne
ylabel ne_z
grid on

%% Fit: 'ne(r)-Bessel'.
[xData, yData] = prepareCurveData( r_ne, ne_r );

% Set up fittype and options.
ft = fittype( 'besselj(0,k2*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.MaxIter = 1000;
opts.Robust = 'Bisquare';
opts.StartPoint = 0.04;

% Fit model to data.
[fitresult{4}, gof(4)] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'ne(r)-Bessel' );
h = plot( fitresult{4}, xData, yData );
legend( h, 'ne_r vs. r_ne', 'ne(r)-Bessel', 'Location', 'NorthEast' );
% Label axes
xlabel r_ne
ylabel ne_r
grid on

%% Fit: 'ne(r)-polyfit'.
[xData, yData] = prepareCurveData( r_ne, ne_r );

% Set up fittype and options.
ft = fittype( 'poly4' );

% Fit model to data.
[fitresult{5}, gof(5)] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', 'ne(r)-polyfit' );
h = plot( fitresult{5}, xData, yData );
legend( h, 'ne_r vs. r_ne', 'ne(r)-polyfit', 'Location', 'NorthEast' );
% Label axes
xlabel r_ne
ylabel ne_r
grid on

%% Fit: 'Te(r)-polyfit'.
[xData, yData] = prepareCurveData( r_Te, Te_r );

% Set up fittype and options.
ft = fittype( 'poly4' );

% Fit model to data.
[fitresult{6}, gof(6)] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', 'Te(r)-polyfit' );
h = plot( fitresult{6}, xData, yData );
legend( h, 'Te_r vs. r_Te', 'Te(r)-polyfit', 'Location', 'NorthEast' );
% Label axes
xlabel r_Te
ylabel Te_r
grid on


