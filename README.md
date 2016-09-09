# SPECFEM3D_Cartesian for ASKI

### Extension package for [ASKI](https://github.com/seismology-RUB/ASKI): use [SPECFEM3D_Cartesian](https://github.com/geodynamics/specfem3d) for solving the seismic forward problem

SPECFEM3D_Cartesian for ASKI, as well as ASKI and some of its components, 
documentation and examples are available under terms of the 
[GNU General Public License](https://github.com/seismology-RUB/ASKI/blob/master/LICENSE)
(version 2 or higher) via [github](https://github.com/seismology-RUB). 
Please find contact addresses [there](https://github.com/seismology-RUB), or visit 
http://www.rub.de/aski in case you want to get in touch with the authors. If you 
encounter any problems installing or using the software, it will be helpful to 
open (or add to) an "issues" topic at the respective repository on gitHub (e.g.
[here for the ASKI main package](https://github.com/seismology-RUB/ASKI/issues))

The main author is Florian Schumacher, Ruhr-University Bochum, Germany. 


## Usage, Documentation

Read the [manual](doc/SPECFEM3D_Cartesian_for_ASKI_manual.pdf) for information on 
how to set parameters correctly and how to use SPECFEM3D_Cartesian as a forward 
solver for ASKI. 


## Requirements

1. You require an installation of the [ASKI main package](https://github.com/seismology-RUB/ASKI):
   ```
   git clone --depth 1 --branch master https://github.com/seismology-RUB/ASKI
   ```
   
   The directory created by the git clone command will be referred to below as `ASKI/`
2. You need a functioning installation of the SPECFEM3D_Cartesian code, including 
   modifications for usage with ASKI:
   * You can either use the basic extract from the SPECFEM3D_Cartesian master branch
     (by November 2015) that comes with this package in subdirectory [specfem3d/](specfem3d/)
   * or use your running installation of SPECFEM3D_Cartesian and extend it for usage
     with ASKI (see "Extend regular SPECFEM3D" below).
     
   Also refer to the [manual](doc/SPECFEM3D_Cartesian_for_ASKI_manual.pdf), sections 1.3, 1.4.
   
   *In both cases, you still must install this package* (item "Installation" below).
       
3. You need basic experience in using the regular SPECFEM3D_Cartesian software!


## Installation

You should clone the latest version of the master branch of the 
[gitHub repository](https://github.com/seismology-RUB/SPECFEM3D_Cartesian_for_ASKI) 
to *the same* directory where you have cloned the ASKI main package to (in the 
[ASKI documentation](https://github.com/seismology-RUB/ASKI/blob/master/doc/ASKI_manual.pdf)
exemplarily called `/your/programs/`). That is, make sure that the git clone command
```
git clone --depth 1 --branch master https://github.com/seismology-RUB/SPECFEM3D_Cartesian_for_ASKI
```

creates the directory `/your/programs/SPECFEM3D_Cartesian_for_ASKI` (also referred to 
below simply as `SPECFEM3D_Cartesian_for_ASKI/`) and ASKI was installed to directory
`/your/programs/ASKI` (also referred to below simply as `ASKI/`) .

You need to compile few more ASKI binaries:
* in [SPECFEM3D_Cartesian_for_ASKI/Makefile](Makefile), set `COMPILER` appropriately, 
  adjust `FFLAGS` if required and set the variables `BLAS`, `LAPACK`, just as you did 
  for installing the ASKI main package
* Execute command
  ```
  make all
  ```
  
  from path `SPECFEM3D_Cartesian_for_ASKI/`.
  
After that, `ASKI/bin/` should contain the new binaries.


## Extend regular SPECFEM3D

If you have a regular SPECFEM3D_Cartesian installation which has the required 
functionality, you may modify your installation for usage with ASKI in the following way.
This procedure was tested for SPECFEM3D_Cartesian (master branch by 7 Nov 2015), extended by 
two important modifications which were commited to the devel branch on 3 september 2015, 
or are about to be commited by the developers team:

1. [src/specfem3D/setup_sources_receivers.f90](specfem3d/src/specfem3D/setup_sources_receivers.f90), 
   subroutine `setup_sources()`, line 180 :<br>
   removing `USE_FORCE_POINT_SOURCE .or.` from the if-clause, i.e. execute
   (re)definition of `t0` only in case of `USE_RICKER_TIME_FUNCTION == .true.`
2. [src/specfem3D/compute_add_sources_viscoelastic.f90](specfem3d/src/specfem3D/compute_add_sources_viscoelastic.f90):<br>
   always call function `comp_source_time_function_gauss()` with half duration `hdur_gaussian(isource)`
   instead of fixed value of `5.d0*DT`
   
If your regular SPECFEM3D_Cartesian installation has this functionality, you can 
extend it for ASKI by the following 12 steps:

1. install SPECFEM3D_Cartesian on your system and make it run, gain 
   experience in using it (below, the installation path is refered to as 
   `specfem3d/`)
2. copy file [SPECFEM3D_Cartesian_for_ASKI/specfem3D_for_ASKI.f90](specfem3D_for_ASKI.f90) to 
   `specfem3d/src/specfem3D/`
3. replace file `specfem3d/src/generate_databases/model_external_values.f90` by 
   [SPECFEM3D_Cartesian_for_ASKI/model_external_values.f90](model_external_values.f90)
4. append content of file [SPECFEM3D_Cartesian_for_ASKI/parallel_ASKI.f90](parallel_ASKI.f90) 
   to file `specfem3d/src/shared/parallel.f90`
5. append content of file [SPECFEM3D_Cartesian_for_ASKI/specfem3D_par_ASKI.f90](specfem3D_par_ASKI.f90)
   to file `specfem3d/src/specfem3D/specfem3D_par.f90`
6. in `specfem3d/src/specfem3D/rules.mk`:<br>
   add the following line into the definition of `specfem3D_OBJECTS` (e.g. before line with `$(EMPTY_MACRO)`)
   ```
   $O/specfem3D_for_ASKI.spec.o \
   ```
   
   (be aware that the above line *must* start with an actual TAB character in order to conform to the GNU-make syntax)
7. in `specfem3d/src/specfem3D/prepare_timerun.F90` in subroutine `prepare_timerun`:<br>
   add the following line at the end of the subroutine, before the statistics output is written to main output file by rank 0:
   ```
   call prepare_timerun_ASKI()
   ```
   
8. in `specfem3d/src/specfem3D/iterate_time.F90` in subroutine `iterate_time`: <br>
   add the following line just before the "enddo" of the time loop
   ```
   call write_ASKI_output()
   ```
   
9. in `specfem3d/src/specfem3D/finalize_simulation.f90` in subroutine `finalize_simulation`: <br>
   add the following line just before the main output file is closed at the end of the subroutine
   ```
   call save_ASKI_output()
   ```
   
10. Set `USE_SOURCES_RECVS_Z = .true.` in `specfem3d/setup/constants.h` (or wherever 
    your file constants.h is located).
11. recompile all SPECFEM3D binaries by executing
    ```
    make
    ```
    
    in directory `specfem3d/`
12. in order to produce ASKI output in SPECFEM3D simulations, copy file 
    [SPECFEM3D_Cartesian_for_ASKI/Par_file_ASKI](Par_file_ASKI) to your respective `DATA/` path
    (which is e.g. `specfem3d/EXAMPLES/my_example/DATA/`, or `specfem3d/DATA/`). This 
    file must be adjusted for any specific simulation (just as all other parameter files), 
    refer to the [manual](doc/SPECFEM3D_Cartesian_for_ASKI_manual.pdf) on how to use it.

Additionally, you may refer to the [manual](doc/SPECFEM3D_Cartesian_for_ASKI_manual.pdf)
section 1.4 for details on extending a regular SPECFEM3D code copy for use with ASKI.

