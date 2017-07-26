clear all

%% This code calculates the relative error in % on the final calculation of
% the TRADITIONAL Udis, depending on the error associated with each of the
% following variables Udis depends on:
%
%           - PDEN_a; average density of shelf water
%           - PDEN_c; average density within the current
%           - H; frontal depth
%           - L; current width

%% Assigns the errors to each variables

DPDEN_a = 3e-2;
DPDEN_c = 2e-2;
DH = 2;     % 0.25
DL = 500;  % 500

% Assign general order of magnitude to each variables:
PDEN_a = 1023.52;
PDEN_c = 1022.78;
H = 20;
L = 58000;

% Calculates the realtive error (DUdis/Udis)

DPEN = PDEN_a-PDEN_c;
rel_error = DPDEN_c.^2 * (-1/DPEN).^2 +...
    DPDEN_a.^2 * (PDEN_c/(PDEN_a*DPEN)).^2 +...
    DH.^2 * (1/H).^2 +...
    DL.^2 * (-1/L).^2;

rel_error = sqrt(rel_error)
    







