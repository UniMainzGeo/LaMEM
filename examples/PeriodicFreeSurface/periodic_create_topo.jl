using GeophysicalModelGenerator

# Creates the initial topography for periodic free surface equilibration test

topo_file="topo.bin"

sigma = 0.1
mu    = 1.5

x_vals = collect(0.0:0.025:2.0)  # x-coordinates
y_vals = collect(0.0:0.025:1.0)  # y-coordinates

X_t, Y_t, Z_t = xyz_grid(x_vals, y_vals, 0)

Topo_z = zeros(length(x_vals), length(y_vals), 1)
    
for i in eachindex(x_vals)
	
	x = x_vals[i]
		
	z_surf = 0.5 + 0.1* exp(-(x-mu)^2/2/sigma^2)/sigma/sqrt(2*pi)	
	
    	Topo_z[i, :, 1] .= z_surf
    	
    	if(x == 2.0)
    		println(x, z_surf)
    	end
    	
    	if(x == 0.0)
    		println(x, z_surf)
    	end
    	
end

Topo = CartData(X_t, Y_t, Z_t, (Topography=Topo_z,))

save_LaMEM_topography(Topo, topo_file)
