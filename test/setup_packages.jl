# this downloads the required packages

# Add PETSc
using Pkg
Pkg.add(name="PETSc_jll", version="3.19.6")

# Copy the relevant directories over
using PETSc_jll

# copy the contents of all directories in a single one
for path in PETSc_jll.PATH_list
    cur_dir = path[1:end-3]   

    # copy mpi directories - we somehow have to do that one by one
    dirs = ["bin","lib","include","share"]
    for d in dirs
        if isdir(joinpath(cur_dir,d))
            run(`sudo -E cp -r $cur_dir/$d /workspace/destdir/`)   
        end 
    end
end

"""
    copy all files 
"""
function cp_files(srcdir, destdir; force=true)
    for f in readdir(srcdir)
        if isfile(joinpath(srcdir,f))
            src = joinpath(srcdir,f)
            dst = joinpath(destdir,f)
            #cp(src, dst, force=force)
            run(`sudo -E cp -r $src $dst`)

        end
    end
    return nothing
end

# And all required dynamic libraries (except petsc)
for srcdir in PETSc_jll.LIBPATH_list
    if !contains(srcdir,"petsc")
        dest_dir = "/workspace/destdir/lib/"
        cp_files(srcdir, dest_dir)
    end
end

# copy PETSc directories
#run(`sudo -E cp -rf $petsc_dir/lib /workspace/destdir`)

# print
run(`ls /workspace/destdir/lib`);



