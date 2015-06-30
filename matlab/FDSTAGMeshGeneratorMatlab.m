function [X,Y,Z,xcoor,ycoor,zcoor, Xpart, Ypart, Zpart] = FDSTAGMeshGeneratorMatlab(npartx,nparty,npartz,fname, RandomNoise)
% This function creates a variable mesh for LaMEM setups based on a
% processor configuration
%    fname - name of the file with the processor configuration;
%               in LaMEM this is saved with -SavePartitioning 1

[Nprocx,Nprocy,Nprocz,xc,yc,zc,xcoor,ycoor,zcoor] = GetProcessorPartitioning(fname);

dx = xcoor(2:end)-xcoor(1:end-1);
dy = ycoor(2:end)-ycoor(1:end-1);
dz = zcoor(2:end)-zcoor(1:end-1);

x = [];
y = [];
z = [];

% create uniform distribution of markers/cell
for i=1:length(dx)
    for j=1:npartx
        if (j==1) & (i==1)
            a = dx(i)/npartx*0.5;
        elseif j==1 & (i~=1)
            a = x(end) + dx(i-1)/npartx*0.5 + dx(i)/npartx*0.5;
        else
            a = x(end)+dx(i)/npartx;
        end
        x(end+1) = a;
    end
end

for i=1:length(dy)
    for j=1:nparty
        if (j==1) & (i==1)
            a = dy(i)/nparty*0.5;
        elseif j==1 & (i~=1)
            a = y(end) + dy(i-1)/nparty*0.5 + dy(i)/nparty*0.5;
        else
            a = y(end)+dy(i)/nparty;
        end
        y(end+1) = a;
    end
end

for i=1:length(dz)
    for j=1:npartz
        if (j==1) & (i==1)
            a = dz(i)/npartz*0.5;
        elseif j==1 & (i~=1)
            a = z(end) + dz(i-1)/npartz*0.5 + dz(i)/npartz*0.5;
        else
            a = z(end)+dz(i)/npartz;
        end
        z(end+1) = a;
    end
end

% shift domain relative to first coord
x = x + xcoor(1);
y = y + ycoor(1);
z = z + zcoor(1);

% create mesh grid
[X,Y,Z] =   meshgrid(x,y,z);


% Add random noise to the particles
dx_vec  =   repmat(dx(:)',[npartx 1]); dx_vec = dx_vec(:)/npartx;
dy_vec  =   repmat(dy(:)',[nparty 1]); dy_vec = dy_vec(:)/nparty;
dz_vec  =   repmat(dz(:)',[npartz 1]); dz_vec = dz_vec(:)/npartz;

if RandomNoise
    [dXNoise,dYNoise,dZNoise] =   meshgrid(dx_vec,dy_vec,dz_vec);
    
    dXNoise = dXNoise.*(rand(size(dXNoise))-0.5);
    dYNoise = dYNoise.*(rand(size(dYNoise))-0.5);
    dZNoise = dZNoise.*(rand(size(dZNoise))-0.5);
    
    Xpart       =   X + dXNoise;
    Ypart       =   Y + dYNoise;
    Zpart       =   Z + dZNoise;
else
    Xpart = X;
    Ypart = Y;
    Zpart = Z;
end

end

% --------------------------------------
function [Nprocx,Nprocy,Nprocz,xc,yc,zc,xcoor,ycoor,zcoor] = GetProcessorPartitioning(test)
% Read Processor Partitioning
fid=PetscOpenFile(test);

Nprocx=read(fid,1,'int32');
Nprocy=read(fid,1,'int32');
Nprocz=read(fid,1,'int32');

nnodx=read(fid,1,'int32');
nnody=read(fid,1,'int32');
nnodz=read(fid,1,'int32');

ix=read(fid,Nprocx+1,'int32');
iy=read(fid,Nprocy+1,'int32');
iz=read(fid,Nprocz+1,'int32');

CharLength=read(fid,1,'double');

xcoor=read(fid,nnodx,'double');
ycoor=read(fid,nnody,'double');
zcoor=read(fid,nnodz,'double');

close(fid);

% Dimensionalize
xcoor = xcoor*CharLength;
ycoor = ycoor*CharLength;
zcoor = zcoor*CharLength;

ix = ix+1;
iy = iy+1;
iz = iz+1;

xc = xcoor(ix);
yc = ycoor(iy);
zc = zcoor(iz);

end