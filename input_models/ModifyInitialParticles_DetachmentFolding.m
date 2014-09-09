%ModifyInitialParticles_DetachmentFolding
%
% (1) Reads in the initial particles (for each processor separately)
% (2) Modifies the particle phases
% (3) Writes the particle information back to file (again for each processor)
%
%

clear

Directory = 'InitialParticles_DetFolding';

% Read in the initial particles on rank 0 (to find the total # of procs)
[num_proc] = ReadInitialParticles(0,Directory);


figure(1),clf
style = {'r.','b.','k.','y.','g.','m.','c.','bo'}
for rank=0:num_proc-1     % loop over all processors
    
    % Read particles
    [num_proc, Part, Partic, info] = ReadInitialParticles(rank,Directory);
    
    % Modify the particles (based on [x,y,z] coordinates of the particles)
    ind = find(Part.x>0 & Part.x<0.2 & Part.z<0.5);
    Part.phase(ind) = 0;
    
    
    % Plot the particles
    num_phases = max(Part.phase);
    for iphase=0:num_phases
        ind= find(Part.phase==iphase);
        style_number = iphase+1;
        
        plot3(Part.x(ind), Part.y(ind), Part.z(ind), style{style_number});
        hold on
        
    end
    axis equal
    
    
    % Save the info of that rank into a LaMEM-readable file
     WriteInitialParticles(rank,Part, Partic, info);
    
end

