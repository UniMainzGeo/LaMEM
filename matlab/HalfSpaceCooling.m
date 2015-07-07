function [T] = HalfSpaceCooling(varargin)
% HalfSpaceCooling
%
% Computes a half-space cooling profile (for an oceanic plate)
%
% [T] = HalfSpaceCooling(z,Age_Myrs, T_mantle, [T_surface] )
%
%       Input parameters:
%
%           z           -  Vector or array with distance to surface in km [positive
%                            is deeper]
%
%           Age_Myrs    -   Thermal age of oceanic plate in Myrs
%           T_mantle    -   Temperature of the mantle
%           T_surface   -   [Optional] Temperature of surface
%
%       Output parameters:
%           
%           T           -   Temperature


z_km        =   varargin{1};
Age_Myrs    =   varargin{2};
T_Mantle    =   varargin{3};
if nargin==4
    T_surface = varargin{4};
else
    T_surface = 0;
end

if nargin>4
    error('Unknown # of input parameters')
end


z_m         =   z_km*1e3;
kappa       =   1e-6;
SecYear     =   3600*24*365.25;


% T_adiabat   =   0.5*z_km;
% T           =   T_adiabat + (T_Mantle-T_adiabat(1)).*erf(z_m./sqrt(kappa*Age_Myrs*SecYear*1e6));

T           =   T_surface + (T_Mantle-T_surface).*erf(z_m./sqrt(kappa*Age_Myrs*SecYear*1e6));


