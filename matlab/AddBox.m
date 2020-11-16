function [Phase,Temp] = AddBox(Phase,Temp,X,Y,Z,BlockBounds, PhaseBlock, varargin)
% Adds a Block to the setup
%
%
%   Usage:
%
%   [Phase,Temp] = AddBox(Phase,Temp,X,Y,Z,BlockBounds, PhaseBlock, varargin)
%
%       Required input parameters:
%           Phase, Temp     :   3D arrays with Phase & Temp [C] info
%           X,Y,Z           :   3D arrays with coordinates
%           BlockBounds     :   [Xleft, Xright, Yfront, Yback, Zbot, Ztop] coordinates of block  
%           PhaseBlock      :   Phase of block
%
%    	Optional input parameter pairs:
%
%           'DipAngle'      :     The angle by which the box dips (w.r.t. middle of upper left side)
%           'StrikeAngle'   :     The strike angle of the box     (w.r.t. middle of upper left side)
%           'RotationPoint' :     Point around which we rotate
%
%           'TempType'      :       {'Constant','Linear','Halfspace','SpreadingRate'}
%           
%                Depending what the 'TempType' options are, we need/can define more parameters:
%
%                'None'     :       Nothing [default]
%
%                'Constant' :       Constant temperature with options:
%                   'cstTemp'       - temperature [Celcius]     
%
%                'Linear' :         Linear temperature profile with options:
%                   'topTemp'       - temperature @ top    [Celcius]
%                   'botTemp'       - temperature @ bottom [Celcius]
%
%                'Halfspace' :      Halfspace cooling temperature profile with options:
%                   'topTemp'       - temperature @ top    [Celcius]
%                   'botTemp'       - temperature @ bottom [Celcius]
%                   'thermalAge'    - thermal age of halfspace cooling profile [Myrs]
%
%                'SpreadingRate' :  Mid Oceanic Ridge with spreading rate
%                   'MORside'       - side @ which the ridge is {'Left','Right','Front','Back'}
%                   'SpreadingVel'  - Spreading velocity away from ridge [cm/yr]
%                   'AgeRidge'      - Thermal age of the ridge [0 by default]
%
%

% Process input parameters & set default parameters (see function below)
p   =   [];
[p] =   ParseInput(p,varargin,'TempType',       'None'	);   
[p] =   ParseInput(p,varargin,'cstTemp',        1000  	);   % Temperature in box in case it has a constant T
[p] =   ParseInput(p,varargin,'botTemp',        1350 	);   % Temperature @ bottom of box
[p] =   ParseInput(p,varargin,'topTemp',        20   	);   % Temperature @ top of plate
[p] =   ParseInput(p,varargin,'thermalAge',     50  	);   % Thermal age of plate if it has a halfspace cooling profile
[p] =   ParseInput(p,varargin,'MORside',        'Left'	);   % Side where the ridge is
[p] =   ParseInput(p,varargin,'SpreadingVel',   5       );   % Spreading velocity away from ridge [cm/yr]
[p] =   ParseInput(p,varargin,'AgeRidge',       0       );   % Thermal age of the ridge [0 by default]
[p] =   ParseInput(p,varargin,'DipAngle',       0       );   % the dip angle of the whole box in degrees
[p] =   ParseInput(p,varargin,'RotationPoint', [0 0 0]  );   % the point that acts as origin around which we rotate 
[p] =   ParseInput(p,varargin,'StrikeAngle',    0       );   % the strike angle of the whole box in degrees


% Find block indices 
Xleft       =   BlockBounds(1);
Xright      =   BlockBounds(2);
Yfront      =   BlockBounds(3);
Yback       =   BlockBounds(4);
Zbot        =   BlockBounds(5);
Ztop        =   BlockBounds(6);

Origin      =   p.Results.RotationPoint;

% Perform rotation if required:
roty     =  @(t) [cosd(-t) 0 sind(-t) ; 0 1 0 ; -sind(-t) 0  cosd(-t)] ;        % rotation around y-axis
rotz     =  @(t) [cosd(t) -sind(t) 0 ; sind(t) cosd(t) 0 ; 0 0 1] ;         % rotation around z-axis
CoordVec =  [X(:)'-Origin(1); Y(:)'-Origin(2); Z(:)'-Origin(3)];   % vector with coordinates
CoordRot =  rotz(p.Results.StrikeAngle)*CoordVec;                           % Rotate around z-axis
CoordRot =  roty(p.Results.DipAngle)*CoordRot;                              % Rotate around y-axis

% transform back to 3D
Xrot = zeros(size(X)); Xrot(find(Xrot==Xrot))= CoordRot(1,:);
Yrot = zeros(size(Y)); Yrot(find(Yrot==Yrot))= CoordRot(2,:);
Zrot = zeros(size(Z)); Zrot(find(Zrot==Zrot))= CoordRot(3,:);

X = Xrot; Y = Yrot; Z = Zrot;

indB        =   find( X>=(Xleft  - Origin(1)) & X<=(Xright - Origin(1)) & ...
                      Y>=(Yfront - Origin(2)) & Y<=(Yback  - Origin(2)) & ...
                      Z>=(Zbot   - Origin(3)) & Z<=(Ztop   - Origin(3)) );

% Set phase of block
Phase(indB) =   PhaseBlock;


% Set temperature structure
switch lower(p.Results.TempType)
    case 'none'
        % do nothing
        
    case 'constant'
        Temp(indB) = p.Results.cstTemp;
         
    case 'linear'
        Ztop = Ztop   - Origin(3);
        Zbot = Zbot   - Origin(3);
        % Linear T profile between Top & bottom
        Tlin = (Z(indB)-Ztop)/(Zbot-Ztop)*(p.Results.botTemp-p.Results.topTemp) + p.Results.topTemp;      
        
        Temp(indB) = Tlin;
    
    case 'halfspace'
        Ztop = Ztop   - Origin(3);
         
        % Set halfspace cooling profile
        kappa       =   1e-6;
        SecYear     =   3600*24*365;
        ThermalAge  =   p.Results.thermalAge*1e6*SecYear;
        T_surface   =   p.Results.topTemp;
        T_mantle    =   p.Results.botTemp;
        
        Temp(indB)  =   (T_surface-T_mantle)*erfc((abs(Z(indB)-Ztop)*1e3)./(2*sqrt(kappa*ThermalAge))) + T_mantle;
        
        
    case 'spreadingrate'
        % Spreading rate thermal structure, which has a MOR on one side
        % (right, left, front, back)
        
        % Shift, in case of rotated box
        Ztop    = Ztop   - Origin(3);
        Xleft   = Xleft  - Origin(1);
        Yfront  = Yfront - Origin(2);
        
        % compute the distance to the ridge
        switch lower(p.Results.MORside)
            case 'left'
                Distance = X(indB)-Xleft;    % in km
            case 'right'
                Distance = Xright-X(indB);   % in km
            case 'front'
                Distance = Y(indB)-Yfront;   % in km
            case 'back'
                Distance = Yback-Y(indB);    % in km
            otherwise 
                error('unknown side')
        end
      
        % Translate distance into thermal age
        ThermalAge  =   (Distance*1e3*1e2)/p.Results.SpreadingVel + p.Results.AgeRidge*1e6;   % Thermal age in years
        SecYear     =   3600*24*365;
        ThermalAge  =   ThermalAge*SecYear;                                                   % in s
        
        T_surface   =   p.Results.topTemp;
        T_mantle    =   p.Results.botTemp;
        kappa       =   1e-6;
        
        % Set halfspace cooling temperature
        Temp(indB)  =   (T_surface-T_mantle).*erfc(abs((Z(indB)-Ztop)*1e3)./(2*sqrt(kappa*ThermalAge))) + T_mantle;
       
        
        
    otherwise
        error('Unknown TempType')
end




% -------------------------------------------------------------------------
function [p] = ParseInput(p,varargin,VariableName,   Value)
% scans through the input parameters


if ~isfield(p,'Results')
    p.Results = [];
end
            
Found = false;
for i=1:2:length(varargin) % scan parameters & find VariableName
    if strcmp(varargin{i},VariableName)
        % found in argument list; set corresponding value
        Found       =   true;
        p.Results   =   setfield(p.Results, varargin{i},varargin{i+1});
    end
end

if ~Found
    % set default
    p.Results   =   setfield(p.Results,VariableName,Value);
end



