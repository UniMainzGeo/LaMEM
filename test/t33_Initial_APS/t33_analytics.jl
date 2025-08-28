using Statistics

function compare_APS(dir, param_file)
    cur_dir = pwd();
    cd(dir)
    outfile0 = "test_markers_no_APS"
    outfile1 = "test_markers_APS_0p5"

    # Tests must be run on 1 core as markers were created with 1 core
    out0 = run_lamem_local_test(param_file, 1, "-mark_load_file ./mdb -out_pvd 1 -out_file_name $outfile0 -out_ptr 1 -out_ptr_APS 1", opt=true)
    data0, t0 = read_LaMEM_timestep(outfile0, 2, pwd(), fields=("plast_strain [ ]",), surf=false, last=true) 
    
    out1 = run_lamem_local_test(param_file, 1, "-mark_load_file ./mdb_APS_0p5 -out_pvd 1 -out_file_name $outfile1 -out_ptr 1 -out_ptr_APS 1", opt=true) 
    data1, t1 = read_LaMEM_timestep(outfile1, 2, pwd(), fields=("plast_strain [ ]",), surf=false, last=true) 
    
    # get the mean value of plast_strain from each model
    APS0 = data0.fields.plast_strain
    APS1 = data1.fields.plast_strain
    cd(cur_dir)
    return mean(APS0), mean(APS1)
end
