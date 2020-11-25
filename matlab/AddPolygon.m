function [Phase,Temp] = AddPolygon(Phase,Temp,X,Y,Z, PolyCoords, PhasePoly, varargin)
% Adds a Block to the setup
%
%
%   Usage:
%
%   [Phase,Temp] = AddPolygon(Phase,Temp,X,Y,Z, PolyCoords, PhasePoly, varargin)
%
%       Required input parameters:
%           Phase, Temp     :   3D arrays with Phase & Temp [C] info
%           X,Y,Z           :   3D arrays with coordinates
%           PolyCoords      :   [Xcoords, ZCoords] coordinates of polygon 
%           PhaseBlock      :   Phase of block
%
%    	Optional input parameter pairs:
%
%           'TempType'      :       {'Constant'}
%           
%                Depending what the 'TempType' options are, we need/can define more parameters:
%
%                'None'     :       Nothing [default]
%
%                'Constant' :       Constant temperature with options:
%                   'cstTemp'       - temperature [Celcius]     
%


% Process input parameters & set default parameters (see function below)
p   =   [];
[p] =   ParseInput(p,varargin,'TempType',           'None'	);   
[p] =   ParseInput(p,varargin,'cstTemp',            1000  	);   % Temperature in box in case it has a constant T


% Polygon coordinates
PolyX          =   PolyCoords(:,1);
PolyZ          =   PolyCoords(:,2);

% find indices inside polygon
indPoly        =   inpolygon(X,Z,PolyX,PolyZ);

% Set phase of polygon
Phase(indPoly) =   PhasePoly;


% Set temperature structure
switch lower(p.Results.TempType)
    case 'none'
        % do nothing
        
    case 'constant'
        Temp(indPoly) = p.Results.cstTemp;
        
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



