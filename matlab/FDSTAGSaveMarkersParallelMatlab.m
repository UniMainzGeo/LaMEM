function FDSTAGSaveMarkersParallelMatlab(varargin)
% This function saves the LaMEM model setup to a parallel marker grid
%   
% FDSTAGSaveMarkersParallelMatlab(A, fname, [Is64BIT])  
%
%  A - structure that contains:
%            Xpart   - 3D matrix containing the x coordinates of markers
%            Ypart   - 3D matrix containing the y coordinates of markers
%            Zpart   - 3D matrix containing the z coordinates of markers
%            x       - 1D vector containing the x coordinates of markers
%            y       - 1D vector containing the y coordinates of markers
%            z       - 1D vector containing the z coordinates of markers
%            Phase   - phase information of the markers
%            Temp    - Temperature information of the markers
%
%   fname   - name of the file with the processor configuration
%             (in LaMEM this is created with -mode save_grid)
%            if your setup is generated on 1 core, pass fname=[]
%
%   Is64BIT - [optional] logical(1) if you are reading a 64 bit file

% ----------- Function begin ----------- %

A       = varargin{1};
fname   = varargin{2};
if nargin==3
    Is64BIT = varargin{3};
else
    Is64BIT = logical(0);
end

if ~isfolder('markers')
    mkdir('markers')
end

% No. of properties the markers carry: x,y,z-coord, phase, T
num_prop = 5;

if isempty(fname)
    % in case we run this on 1 processor only
    Nprocx  = 1;
    Nprocy  = 1;
    Nprocz  = 1;
    xc      = A.x;
    yc      = A.y;
    zc      = A.z;
    
else
    
    % Read Processor Partitioning
    [P] = GetProcessorPartitioning(fname, Is64BIT);
    
    % get number of processors and processor coordnate bounds
    Nprocx = P.Nprocx;
    Nprocy = P.Nprocy;
    Nprocz = P.Nprocz;
    xc     = P.xc;
    yc     = P.yc;
    zc     = P.zc;
end

Nproc                      = Nprocx*Nprocy*Nprocz;
[num, num_i, num_j, num_k] = get_numscheme(Nprocx, Nprocy, Nprocz);

% Particle coordinates (should be permuted such that it has the same size
% as A.phase)
X = A.Xpart;
Y = A.Ypart;
Z = A.Zpart;

% Get particles of respective procs
[xi,ix_start,ix_end] = get_ind(A.x,xc,Nprocx);
[yi,iy_start,iy_end] = get_ind(A.y,yc,Nprocy);
[zi,iz_start,iz_end] = get_ind(A.z,zc,Nprocz);

x_start(num)= ix_start(num_i);
x_end(num)  = ix_end(num_i);
y_start(num)= iy_start(num_j);
y_end(num)  = iy_end(num_j);
z_start(num)= iz_start(num_k);
z_end(num)  = iz_end(num_k);

% Loop over all processors partition

for num=1:Nproc

    part_x   = X(x_start(num):x_end(num),y_start(num):y_end(num),z_start(num):z_end(num));
    part_y   = Y(x_start(num):x_end(num),y_start(num):y_end(num),z_start(num):z_end(num));
    part_z   = Z(x_start(num):x_end(num),y_start(num):y_end(num),z_start(num):z_end(num));
    part_phs = A.Phase(x_start(num):x_end(num),y_start(num):y_end(num),z_start(num):z_end(num));
    part_T   = A.Temp(x_start(num):x_end(num),y_start(num):y_end(num),z_start(num):z_end(num));

    % No. of particles per processor
    num_particles = size(part_x,1)* size(part_x,2) * size(part_x,3);

    % Information vector per processor
    lvec_info(1)  = num_particles;

    lvec_prtcls = zeros(1,num_prop*num_particles);

    lvec_prtcls(1:num_prop:end) = part_x(:);
    lvec_prtcls(2:num_prop:end) = part_y(:);
    lvec_prtcls(3:num_prop:end) = part_z(:);
    lvec_prtcls(4:num_prop:end) = part_phs(:);
    lvec_prtcls(5:num_prop:end) = part_T(:);

    % Output files
    fname = sprintf('./markers/mdb.%1.8d.dat', num-1);
    disp(['Writing file -> ',fname])
    lvec_output    = [lvec_info(:); lvec_prtcls(:)];

    PetscBinaryWrite(fname,lvec_output);

    clear part_x part_y part_z part_phs lvec_info lvec_prtcls id_sort No_id_vec No_id lvec_output
end

end
% ----------- Function end ----------- %

% --------------------------------------
function [n,nix,njy,nkz] = get_numscheme(Nprocx,Nprocy,Nprocz)
num=0;
for k=1:Nprocz
    for j=1:Nprocy
        for i=1:Nprocx
            num=num+1;
            n(num)   = num;
            nix(num)= i;
            njy(num)= j;
            nkz(num)= k;
        end
    end
end

end

% --------------------------------------
function [xi,ix_start,ix_end] = get_ind(x,xc,Nprocx)
    if Nprocx == 1
        xi       = length(x);
        ix_start = 1;
        ix_end   = length(x);
    else
        for k= 1:Nprocx
                if k==1
                    xi(k) = length(x(x>=xc(k)& x<=xc(k+1)));
                else
                    xi(k) = length(x(x>xc(k) & x<=xc(k+1)));
                end
        end
        ix_start = cumsum([0,xi(1:end-1)])+1;
        ix_end   = cumsum(xi(1:end));
    end
end
