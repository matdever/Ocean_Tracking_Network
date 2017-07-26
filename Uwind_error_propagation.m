clear all

%% This code calculates the relative error in % on the final calculation of
% the TRADITIONAL Uwind, depending on the error associated with each of the
% following variables Udis depends on:
%
%           - RHO_air; air density
%           - RHO; depth average ocean density
%           - C10; surface drag coeficient
%           - CDa; depth average drag coeficient
%           - U10; alongshore wind speed at 10m

%% Assigns the errors to each variables

DRHO_air = 0;
DRHO = 1;
DC10 = 0;     % 0.25
DCDa = 0;  % 500
DU10 = .01;

% Assign general order of magnitude to each variables:
RHO_air = 1.2;
RHO = 1024;
C10 = 1.2e-3;
CDa = 1/(8+(1/sqrt(2.5e-3)))^2;
U10 = 30;

% Calculates the realtive error (DUwind/Uwind)

rel_error = DRHO_air.^2 * (1/2/RHO_air).^2 +...
    DRHO.^2 * (-1/2/RHO).^2 + ...
    DC10.^2 * (1/2/C10).^2 + ...
    DCDa.^2 * (-1/2/CDa).^2 + ...
    DU10.^2 * (1/U10).^2;

rel_error = sqrt(rel_error)
    







