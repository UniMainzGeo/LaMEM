# Debugging LaMEM and Microsoft Visual Studio Code

Microsoft Visual Studio Code is a releatively new IDE, which we found quite useful in combination with LaMEM.
Here, we provide some instruction that helped us do this on a Mac. Linux instructions are likely similar (probably easier).


### 1.1.2. Install the debugger and make sure it works on Mac
We assume that you download and install VSC to your machine and installed the C/C++ package as well. 

First, you will have to make sure that you installed the gnu debugger (provided that your version of PETSc was installed with gcc compilers). Important is ofcourse that the versions are compatible. One possibility is to install MacPorts, and install compilers with:
```
$ sudo port install gdb
```
Once installed we have to ensure that it actually works on a mac, as the safety precautions of mac do not allow gdb to control other processes. Therefore, you need to *codesign* gdb by creating a certificate that allows it to hack your executable.

An explanation of how that can be done is given here:

https://sourceware.org/gdb/wiki/BuildingOnDarwin#Giving_gdb_permission_to_control_other_processes%7Cofficial


### 1.1.2. Linux


