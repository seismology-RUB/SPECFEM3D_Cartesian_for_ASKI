#!/usr/bin/python
#
#!/usr/bin/env /usr/bin/python
#
#----------------------------------------------------------------------------
#   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
#
#   This file is part of ASKI version 1.2.
#
#   ASKI version 1.2 is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   ASKI version 1.2 is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------------
#
# import python modules
from os import system as os_system
from os import environ as os_environ
from os import mkdir as os_mkdir
from os import path as os_path
from os import listdir as os_listdir
from sys import exit as sys_exit
from sys import path as sys_path
from time import time as time_time
from time import ctime as time_ctime
#
# get SUN GRID ENGINE environmental variables (if any)
SGE_job_id = os_environ.get('JOB_ID')
runs_on_SGE = SGE_job_id is not None
SGE_o_workdir = os_environ.get('SGE_O_WORKDIR')
SGE_hostname = os_environ.get('HOSTNAME')
SGE_pe_hostfile = os_environ.get('PE_HOSTFILE')
if SGE_pe_hostfile is not None:
    SGE_pe_hostfile_content = open(SGE_pe_hostfile).read()
else:
    SGE_pe_hostfile_content = 'no content, file PE_HOSTFILE is empty'
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# IMPORTANT STUFF TO ADJUST
#
# import own modules
# when running in SUN GRID ENGINE (or maybe also on other HPC queueing systems), the PYTHONPATH
# environment variable must be additionally told, where your own python modules are, before you
# can import them.
# There are two options here:
# 1) manually append the src path of your ASKI main package installation to variable sys_path:
#######sys_path.append("/home/florian/code/ASKI/src")
# 2) copy the files readEventStationFile.py and inputParameter.py to the path from which this
#    script is called. Then append current path "./" to the PYTHONPATH
sys_path.append("./")
# Following either 1) or 2), now the following modules can be imported:
from inputParameter import inputParameter
from readEventStationFile import eventList,stationList
#
# main parameter file of inversion
main_parfile = '/rscratch/minos18/florian/ASKI/ASKI_inversion_cross_borehole/main_parfile_cross_borehole'
#
# define the command which will be called (via system call) for each simulation;
# e.g. can be something like './process_solver_only.sh' along with using a different command in the first iteration (see below)
command_system_call = './process_solver_only.sh'
#
# if in the VERY FIRST iteration (i.e. the first simulation that is done) a DIFFERENT
# command should be issued, indicate so by the following flag, and define the alternative command
use_different_command_in_first_simulation = True
command_system_call_first_simulation = './process.sh'
#
# say, if ASKI output volume, which is used in SPECFEM, should be defined by the inversion grid 
# definition of the current inversion grid (if possible)
define_ASKI_output_volume_by_inversion_grid = True
#
# OUTPUT_FILES_PATH and IN_DATA_FILES_PATH (as defined in SPECFEM3D/src/shared/constants.h, 
# and Par_file, BUT relative to the current path, where this script is run)
OUTPUT_FILES_PATH = 'OUTPUT_FILES/'
LOCAL_PATH = 'OUTPUT_FILES/DATABASES_MPI/'
IN_DATA_FILES_PATH = 'DATA/'
#
# log file name, will be written to OUTPUT_FILES_PATH
logfile = 'run_specfem3dCartesianForASKI_simulations.log'
#
if runs_on_SGE:
    logfile += str(SGE_job_id)
#
# EMAILING LOGFILE SOMEWHERE
#
# set to False if no emails should be sent at all
send_emails = False
# email address, to which log notifications are sent
email_receiver= 'receiver@mail.domain'
email_sender = 'sender@mail.domain'
# if send_emails = True, the script always writes an email after the 1st iteration and at the end of all iterations (or if script exits unintendedly)
# number_of_emails_during_iteration defines the number of additional emails in between (1st and last iteration) while iterating over the simulations
# set number_of_emails_during_iteration = 0 if you do not want to receive any additional emails aside from the two after 1st and last iteration
number_of_emails_during_iteration = 0
#
####################################################################################################################
## define the (order of the) specfem3dForASKI simulations by strings displ_simulations,gt_simulations,measured_data_simulations
## in the following way:
##
## string displ_simulations defines the events for which kernel_displacement output for ASKI is computed, 
## string gt_simulations defines the (components of) stations for which kernel_green_tensor output for ASKI is computed,
## string measurde_data_simuulations defines the events for which synthetically calculated "measured_data" can be computed
##   (w.r.t a perturbed model) in case of pure synthetic test studies, these are regular SPECFEM3D simulations without producing
##   any ASKI output (in general, these simulations are made separately from kernel simulations, as you will use a different
##   model and do not need to compute ASKI output for kernels)
##
## displ_simulations must be of one of the folling forms:
##   displ_simulations = ''
##      no kernel_displacement simulations will be done
##   displ_simulations = 'all'
##      in this case kernel_displacement output is computed for all events defined in FILE_EVENT_LIST
##   displ_simulations = 'Source024,Source002,Source005'
##      for all events in the ','-separated list of eventIDs (here 3 events: Source024, Source002 and Source005), 
##      kernel_displacement output is computed (eventIDs must be present in FILE_EVENT_LIST)
##   displ_simulations = 'all-except:S001,S002,S024'
##      all events are taken into account, except the ones defined by the ','-separated list of eventIDs
##      following the word 'all-except:'  (here all events except the three S001,S002 and S024)
##   IF THERE IS ANY INVALID eventID (i.e. not present in FILE_EVENT_LIST), THIS SCRIPT WILL RAISE AN ERROR!
##
## gt_simulations must be of one of the folling forms:
##   gt_simulations = ''
##      no kernel_green_tensor simulations will be done
##   gt_simulations = 'all'
##      in this case kernel_green_tensor output is computed for all stations defined in FILE_STATION_LIST for 
##      all components defined by gt_components (see below for form of gt_components)
##   gt_simulations = 'all-except:Geophone23,Geophone2'
##      all stations are taken into account (at all components defined by gt_components), except the ones defined by 
##      the ','-separated list of station names following the word 'all-except:' (here all stations excpet Geophone23 and Geophone2)
##   gt_simulations = 'specific'
##      in this case kernel_green_tensor output is computed for specific components of specific stations, BOTH
##      defined by gt_components (see below for form of gt_components)
##   IF THERE IS ANY INVALID station name (i.e. not present in FILE_STATION_LIST), THIS SCRIPT WILL RAISE AN ERROR!
##
## gt_components must be of the following form in case of gt_simulations being 'all' or 'all-except:...' :
##   gt_components = 'CX,CZ'
##     a ','-separated list of valid components (the currently supported components are the global Cartesian 
##     components 'CX','CY','CZ')
##   IF THERE IS ANY INVALID COMPONENT, THIS SCRIPT WILL RAISE AN ERROR!
##
## gt_components must be of the following form in case of gt_simulations being 'specific' :
##   gt_components = 'Geophone23:CX,CY,CZ;Geophone2:CZ;ReceiverB:CY,CX'
##     a ';'-separated list of entries consisting of the station name and ':' followed by a ','-separated list of 
##     valid components. This defines the stations and the station-specific components for which this script
##     computes Green functions. 
##   IF THERE IS ANY INVALID COMPONENT, THIS SCRIPT WILL RAISE AN ERROR!
## 
## measured_data_simulations must be of the same form as displ_simulations (see above)
####################################################################################################################
displ_simulations = 'all'
gt_simulations = 'all'
gt_components = 'CX,CY,CZ'
measured_data_simulations = ''
#
# DATA OUTPUT
#
# specify here, whether or not the SPECFEM3D STATIONS file should be produced by ASKI (based on ASKI's stations file)
create_specfem_stations = True
#
#
# END OF STUFF TO ADJUST
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#
############################################################
# CLASS simulation
############################################################
#
class simulation:

    def __init__(self):

        # number of processes defined in process.sh (needed for function check_process_sh() )
        nproc_list = [line.split('=')[1] for line in open('./process.sh','r') if '=' in line and line.split('=')[0] == 'NPROC']
        if len(nproc_list) != 1:
            self.log("### STOP :  found "+str(len(nproc_list))+" possible values for 'NPROC' in './process.sh': '"+"', '".join(nproc_list)+"'\n")
            raise Exception("not exactly one possible value for NPROC in ./process.sh; see logfile '"+logfile+"'")
        else:
            self.nproc = int(nproc_list[0])

        # open main parfile and check if all required keywords are present
        if not os_path.exists(main_parfile):
            self.log("### STOP : file '"+main_parfile+"' as set for the main parameter file does not exist, "+
                     "please correct the definition of 'main_parfile = ...' at the beginning of this script\n\n")
            raise Exception("main parameter file does not exist; see logfile '"+logfile+"'")
        try:
            self.mparam = inputParameter(main_parfile)
        except:
            self.log("### STOP : could not create inputParameter object for main parameter file '"+
                     main_parfile+"', make sure that the file is of correct form\n\n")
            raise
        # check if all required keys are set
        noKeys = self.mparam.keysNotPresent(['MAIN_PATH_INVERSION','CURRENT_ITERATION_STEP','ITERATION_STEP_PATH',
                                            'PARFILE_ITERATION_STEP','FILE_STATION_LIST','FILE_EVENT_LIST',
                                             'PATH_MEASURED_DATA','MEASURED_DATA_FREQUENCY_STEP'])
        if len(noKeys) > 0:
            self.log("### STOP : the following keywords are required in the main parameter file '"+
                     main_parfile+"':\n"+
                     "### "+',  '.join(noKeys)+"\n\n")
            raise Exception("missing keywords in main paramter file; see logfile '"+logfile+"'")

        # store iteration step path for further use
        self.iter_path = os_path.join(self.mparam.sval('MAIN_PATH_INVERSION'),self.mparam.sval('ITERATION_STEP_PATH')+
                                      '%3.3i/'%self.mparam.ival('CURRENT_ITERATION_STEP'))

        # open iteration step parfile and check if all required keywords are present
        iter_parfile = os_path.join(self.iter_path,self.mparam.sval('PARFILE_ITERATION_STEP'))
        if not os_path.exists(iter_parfile):
            self.log("### STOP : the iteration step parameter file '"+iter_parfile+"' as derived from the main "+
                     "parameter file does not exist, please make sure that the settings of MAIN_PATH_INVERSION, "+
                     "CURRENT_ITERATION_STEP, ITERATION_STEP_PATH and PARFILE_ITERATION_STEP in main parameter file '"+
                     main_parfile+"' is correctly set\n\n")
            raise Exception("iteration step parameter file does not exist; see logfile '"+logfile+"'")
        try:
            self.iparam = inputParameter(iter_parfile)
        except:
            self.log("### STOP : could not create inputParameter object for the iteration step parameter file '"+
                     iter_parfile+"', make sure that the file is of correct form\n\n")
            raise
        # check if all required keys are set
        noKeys = self.iparam.keysNotPresent(['ITERATION_STEP_NUMBER_OF_FREQ','ITERATION_STEP_INDEX_OF_FREQ',
                                             'PATH_KERNEL_DISPLACEMENTS','PATH_KERNEL_GREEN_TENSORS',
                                             'TYPE_INVERSION_GRID','PARFILE_INVERSION_GRID'])
        if len(noKeys) > 0:
            self.log("### STOP : the following keywords are required in iteration step parameter file '"+
                     iter_parfile+"':\n"+
                     "### "+',  '.join(noKeys)+"\n\n")
            raise Exception("missing keywords in iteration step paramter file; see logfile '"+logfile+"'")

        # in case of define_ASKI_output_volume_by_inversion_grid, we also need the content of PARFILE_INVERSION_GRID
        if define_ASKI_output_volume_by_inversion_grid:
            if self.iparam.sval('TYPE_INVERSION_GRID') == 'scartInversionGrid':
                ASKI_type_inversion_grid_char = 'scartInversionGrid'
                parfile_invgrid = os_path.join(self.iter_path,self.iparam.sval('PARFILE_INVERSION_GRID'))
                try:
                    param = inputParameter(parfile_invgrid)
                except:
                    self.log("### STOP : could not create inputParameter object for the 'scartInversionGrid' parameter file '"+
                             parfile_invgrid+"' (for defining ASKI output volume from inversion grid)\n\n")
                    raise
                # check if all required keys are set
                noKeys = param.keysNotPresent(['SCART_INVGRID_CX','SCART_INVGRID_CY','SCART_INVGRID_ZMAX',
                                               'SCART_INVGRID_WX','SCART_INVGRID_WY','SCART_INVGRID_ROT',
                                               'SCART_INVGRID_NREF_BLOCKS','SCART_INVGRID_NLAY',
                                               'SCART_INVGRID_THICKNESS','SCART_INVGRID_NX','SCART_INVGRID_NY'])
                if len(noKeys) > 0:
                    self.log("### STOP : the following keywords are required in 'scartInversionGrid' parameter file '"+
                             self.iparam.sval('PARFILE_INVERSION_GRID')+"' (for defining ASKI output volume from inversion grid):\n"+
                             "### "+',  '.join(noKeys)+"\n\n")
                    raise Exception("missing keywords in 'scartInversionGrid' paramter file (for defining ASKI output volume from inversion grid); see logfile '"+logfile+"'")
                # SCART_INVGRID_WX
                if param.fval('SCART_INVGRID_WX') is not None and param.fval('SCART_INVGRID_WX') > 0.:
                    self.ASKI_wx = param.sval('SCART_INVGRID_WX')
                else:
                    raise Exception("'SCART_INVGRID_WX' = '"+param.sval('SCART_INVGRID_WX')+
                                    "' is not valid in 'scartInversionGrid' parameter file '"+
                                    self.iparam.sval('PARFILE_INVERSION_GRID')+"', must be positive")
                # SCART_INVGRID_WY
                if param.fval('SCART_INVGRID_WY') is not None and param.fval('SCART_INVGRID_WY') > 0.:
                    self.ASKI_wy = param.sval('SCART_INVGRID_WY')
                else:
                    raise Exception("'SCART_INVGRID_WY' = '"+param.sval('SCART_INVGRID_WY')+
                                    "' is not valid in 'scartInversionGrid' parameter file '"+
                                    self.iparam.sval('PARFILE_INVERSION_GRID')+"', must be positive")
                # SCART_INVGRID_CX
                if param.fval('SCART_INVGRID_CX') is not None:
                    self.ASKI_cx = param.sval('SCART_INVGRID_CX')
                else:
                    raise Exception("'SCART_INVGRID_CX' = '"+param.sval('SCART_INVGRID_CX')+
                                    "' is no real number in 'scartInversionGrid' parameter file '"+
                                    self.iparam.sval('PARFILE_INVERSION_GRID')+"'")
                # SCART_INVGRID_CY
                if param.fval('SCART_INVGRID_CY') is not None:
                    self.ASKI_cy = param.sval('SCART_INVGRID_CY')
                else:
                    raise Exception("'SCART_INVGRID_CY' = '"+param.sval('SCART_INVGRID_CY')+
                                    "' is no real number in 'scartInversionGrid' parameter file '"+
                                    self.iparam.sval('PARFILE_INVERSION_GRID')+"'")
                # SCART_INVGRID_ROT
                if param.fval('SCART_INVGRID_ROT') is not None:
                    self.ASKI_rot_Z = param.sval('SCART_INVGRID_ROT')
                else:
                    raise Exception("'SCART_INVGRID_ROT' = '"+param.sval('SCART_INVGRID_ROT')+
                                    "' is no real number in 'scartInversionGrid' parameter file '"+
                                    self.iparam.sval('PARFILE_INVERSION_GRID')+"'")
                # SCART_INVGRID_ZMAX
                if param.fval('SCART_INVGRID_ZMAX') is not None:
                    zmax = param.fval('SCART_INVGRID_ZMAX')
                else:
                    raise Exception("'SCART_INVGRID_ZMAX' = '"+param.sval('SCART_INVGRID_ZMAX')+
                                    "' is no real number in 'scartInversionGrid' parameter file '"+
                                    self.iparam.sval('PARFILE_INVERSION_GRID')+"'")
                # SCART_INVGRID_NREF_BLOCKS
                if param.ival('SCART_INVGRID_NREF_BLOCKS') is not None and param.ival('SCART_INVGRID_NREF_BLOCKS') > 0:
                    nref_blocks = param.ival('SCART_INVGRID_NREF_BLOCKS')
                else:
                    raise Exception("'SCART_INVGRID_NREF_BLOCKS' = '"+param.sval('SCART_INVGRID_NREF_BLOCKS')+
                                    "' is not valid in 'scartInversionGrid' parameter file '"+
                                    self.iparam.sval('PARFILE_INVERSION_GRID')+"', must be positive")
                # SCART_INVGRID_NLAY
                if param.ilist('SCART_INVGRID_NLAY',nref_blocks) is not None:
                    nlay = param.ilist('SCART_INVGRID_NLAY',nref_blocks)
                    if any([ilay <= 0 for ilay in nlay]):
                        raise Exception("all entries of vector 'SCART_INVGRID_NLAY' = '"+param.sval('SCART_INVGRID_NLAY')+
                                        "' must be positive in 'scartInversionGrid' parameter file '"+
                                        self.iparam.sval('PARFILE_INVERSION_GRID'))
                else:
                    raise Exception("'SCART_INVGRID_NLAY' = '"+param.sval('SCART_INVGRID_NLAY')+
                                    "' is not a vector of "+str(nref_blocks)+" integers in 'scartInversionGrid' parameter file '"+
                                    self.iparam.sval('PARFILE_INVERSION_GRID'))
                # SCART_INVGRID_THICKNESS
                if param.flist('SCART_INVGRID_THICKNESS',nref_blocks) is not None:
                    thickness = param.flist('SCART_INVGRID_THICKNESS',nref_blocks)
                    if any([t <= 0. for t in thickness]):
                        raise Exception("all entries of vector 'SCART_INVGRID_THICKNESS' = '"+param.sval('SCART_INVGRID_THICKNESS')+
                                        "' must be positive in 'scartInversionGrid' parameter file '"+
                                        self.iparam.sval('PARFILE_INVERSION_GRID'))
                else:
                    raise Exception("'SCART_INVGRID_THICKNESS' = '"+param.sval('SCART_INVGRID_THICKNESS')+
                                    "' is not a vector of "+str(nref_blocks)+" integers in 'scartInversionGrid' parameter file '"+
                                    self.iparam.sval('PARFILE_INVERSION_GRID'))
                wz = sum([i*t for (i,t) in zip(nlay,thickness)])
                self.ASKI_wz = str(wz)
                self.ASKI_cz = str(zmax-0.5*wz)
                self.ASKI_rot_X = '0.0'
                self.ASKI_rot_Y = '0.0'
                self.ASKI_type_inversion_grid = '2'
            #elif self.iparam.sval('TYPE_INVERSION_GRID') == 'ecartInversionGrid':
            #    pass
            else:
                raise Exception("'TYPE_INVERSION_GRID' = '"+self.iparam.sval('TYPE_INVERSION_GRID')+
                                "' is not supported by the automatic definition of ASKI output volume "+
                                "by the current inversion grid. Please define the ASKI output volume manually and set "+
                                "define_ASKI_output_volume_by_inversion_grid = False on top of this script")

        # read event list and station list
        if displ_simulations != '' or measured_data_simulations != '':
            if not os_path.exists(self.mparam.sval('FILE_EVENT_LIST')):
                self.log("### STOP : the event list file '"+self.mparam.sval('FILE_EVENT_LIST')+"' as set in the "+
                         "main parameter file '"+main_parfile+"' does not exist\n\n")
                raise Exception("event list file does not exist; see logfile '"+logfile+"'")
            try:
                self.evlist = eventList(self.mparam.sval('FILE_EVENT_LIST'),list_type='standard')
            except:
                self.log("### STOP : could not construct event list from file '"+self.mparam.sval('FILE_EVENT_LIST')+
                         "', make sure that the file is of correct form\n\n")
                raise
            if self.evlist.nev == 0:
                self.log("### STOP : event list from file '"+self.mparam.sval('FILE_EVENT_LIST')+"' does not contain"+
                         "any valid events\n\n")
                raise Exception("no events in event list file; see logfile '"+logfile+"'")
            if self.evlist.csys != 'C':
                self.log("### STOP : event list from file '"+self.mparam.sval('FILE_EVENT_LIST')+"' is for coordinate system '"+
                         self.evlist.csys+"', only 'C' supported here\n\n")
                raise Exception("coordinate system of event list file not supported; see logfile '"+logfile+"'")
            if create_specfem_stations:
                if not os_path.exists(self.mparam.sval('FILE_STATION_LIST')):
                    self.log("### STOP : the station list file '"+self.mparam.sval('FILE_STATION_LIST')+"' as set "+
                             "in the main parameter file '"+main_parfile+"' does not exist\n\n")
                    raise Exception("station list file does not exist; see logfile '"+logfile+"'")
                try:
                    self.statlist = stationList(self.mparam.sval('FILE_STATION_LIST'),list_type='standard')
                except:
                    self.log("### STOP : could not construct station list from file '"+self.mparam.sval('FILE_STATION_LIST')+
                             "', make sure that the file is of correct form\n\n")
                    raise
                if self.statlist.nstat == 0:
                    self.log("### STOP : station list from file '"+self.mparam.sval('FILE_STATION_LIST')+"' does not contain"+
                             "any valid stations\n\n")
                    raise Exception("no stations in ASKI station list file; see logfile '"+logfile+"'")
                if self.statlist.csys != 'C':
                    self.log("### STOP : station list from file '"+self.mparam.sval('FILE_EVENT_LIST')+"' is for coordinate system '"+
                             self.statlist.csys+"', only 'C' supported here\n\n")
                    raise Exception("coordinate system of station list file not supported; see logfile '"+logfile+"'")
        if gt_simulations != '' and not hasattr(self,'statlist'):
            if not os_path.exists(self.mparam.sval('FILE_STATION_LIST')):
                self.log("### STOP : the station list file '"+self.mparam.sval('FILE_STATION_LIST')+"' as set "+
                         "in the main parameter file '"+main_parfile+"' does not exist\n\n")
                raise Exception("station list file does not exist; see logfile '"+logfile+"'")
            try:
                self.statlist = stationList(self.mparam.sval('FILE_STATION_LIST'),list_type='standard')
            except:
                self.log("### STOP : could not construct station list from file '"+self.mparam.sval('FILE_STATION_LIST')+
                         "', make sure that the file is of correct form\n\n")
                raise
            if self.statlist.nstat == 0:
                self.log("### STOP : station list from file '"+self.mparam.sval('FILE_STATION_LIST')+"' does not contain"+
                         "any valid stations\n\n")
                raise Exception("no stations in ASKI station list file; see logfile '"+logfile+"'")
            if self.statlist.csys != 'C':
                self.log("### STOP : station list from file '"+self.mparam.sval('FILE_EVENT_LIST')+"' is for coordinate system '"+
                         self.statlist.csys+"', only 'C' supported here\n\n")
                raise Exception("coordinate system of station list file not supported; see logfile '"+logfile+"'")
            
        self.all_tasks = []
        self.gt_comp_files = []
        self.append_valid_tasks()
        for staname,file_content in self.gt_comp_files:
            filename = os_path.join(self.iter_path,self.iparam.sval('PATH_KERNEL_GREEN_TENSORS'),'kernel_gt_'+staname+'.comp')
            try:
                open(filename,'w').write(file_content)
            except:
                self.log("### STOP : could not write Green tensor components file '"+filename+"' for station '"+
                         staname+"'\n\n")
                raise Exception("could not write a Green tensor components file; see logfile '"+logfile+"'")

        self.index_iteration_send_email = [0]
        if type(number_of_emails_during_iteration) is int:
            self.index_iteration_send_email.extend([int( (i+1.)*max(float(len(self.all_tasks)-1)/float(number_of_emails_during_iteration+1),1.) )
                                                    for i in range(min(number_of_emails_during_iteration,len(self.all_tasks)-2))])

        # log initial information
        if runs_on_SGE:
            log_SGE_info = ('running on SUN GRID ENGINE with JOB_ID '+str(SGE_job_id)+'; master host is '+
                            str(SGE_hostname)+'; content of PE_HOSTFILE is \n--- START CONTENT PE_HOSTFILE ---\n'+
                            SGE_pe_hostfile_content+'--- END CONTENT PE_HOSTFILE ---')
        else:
            log_SGE_info = ''
        self.log('################################################################################\n'+
                 "Welcome to these automated SPECFEM3D_Cartesian simulations for ASKI\n"+
                 log_SGE_info+"\n\n"+
                 "main ASKI parameter file: '"+main_parfile+"'\n"+
                 "iteration step %i\n"%self.mparam.ival('CURRENT_ITERATION_STEP')+
                 "iteration step specific parameter file: '"+iter_parfile+"'\n"+
                 "'./process.sh' tells me, we're using "+str(self.nproc)+" procs\n"+
                 "OUTPUT_FILES_PATH = '"+OUTPUT_FILES_PATH+"'\n"+
                 "LOCAL_PATH = '"+LOCAL_PATH+"'\n"+
                 "IN_DATA_FILES_PATH = '"+IN_DATA_FILES_PATH+"'\n"+
                 "\n"+
                 "according to the simulation strings \n"+
                 "  displ_simulations = '"+displ_simulations+"'\n"+
                 "  gt_simulations = '"+gt_simulations+"'\n"+
                 "  gt_components = '"+gt_components+"'\n"+
                 "  measured_data_simulations = '"+measured_data_simulations+"',\n"+
                 "now the following "+str(len(self.all_tasks))+" simulations are done (in this order): \n\n"+
                 "(TYPE ID[_component])\n"+
                 ',  '.join(["(%s %s)"%typ_task for typ_task in self.all_tasks])+"\n"+
                 '\n')
        if len(self.gt_comp_files) > 0:
            self.log("already in advance, all "+str(len(self.gt_comp_files))+" Green tensor component files "+
                     "(containing the Green tensor components for each station) were written to paths '"+
                     os_path.join(self.iter_path,self.iparam.sval('PATH_KERNEL_GREEN_TENSORS'),'kernel_gt_staname.comp')+"'\n"+
                     "\n")
        if send_emails:
            self.log("this log will be sent via email to '"+email_receiver+"' after the following simulations (indices, first index is 1):\n"+
                     ', '.join([str(i+1) for i in self.index_iteration_send_email])+'\n'
                     'and when all simulations are successfully finished (or the script exited unintendedly)\n'+
                     '\n')
        else:
            self.log("this log will not be sent via email anywhere\n\n")

        # info about ASKI output volume
        if define_ASKI_output_volume_by_inversion_grid and (displ_simulations != '' or gt_simulations!= '') :
            self.log("for every simulation of type displ or gt, the ASKI output volume in Par_file_ASKI\n"+
                     "will be defined by the following values:\n"+
                     "  ASKI_type_inversion_grid = "+self.ASKI_type_inversion_grid+" (meaning '"+ASKI_type_inversion_grid_char+"')\n"+
                     "  ASKI_wx = "+self.ASKI_wx+"\n"+
                     "  ASKI_wy = "+self.ASKI_wy+"\n"+
                     "  ASKI_wz = "+self.ASKI_wz+"\n"+
                     "  ASKI_rot_X = "+self.ASKI_rot_X+"\n"+
                     "  ASKI_rot_Y = "+self.ASKI_rot_Y+"\n"+
                     "  ASKI_rot_Z = "+self.ASKI_rot_Z+"\n"+
                     "  ASKI_cx = "+self.ASKI_cx+"\n"+
                     "  ASKI_cy = "+self.ASKI_cy+"\n"+
                     "  ASKI_cz = "+self.ASKI_cz+"\n"+
                     "\n")
        else:
            self.log("the ASKI output volume in Par_file_ASKI will not be defined by this script,\n"+
                     "you should have taken care of that yourself\n\n")

        self.time_start = time_time()
        self.log("starting simulations now at time -- "+time_ctime(self.time_start)+"\n"+
                 "\n")
#
#-----------------------------------------------------------
#
    def append_valid_tasks(self):#displ_simulations,gt_simulations,measured_data_simulations)
        # extract all valid eventIDs and station names for testing, omit numeric keys (event index, station index) present in dictionaries
        if displ_simulations != '' or measured_data_simulations != '':
            valid_evids = sorted([key for key in self.evlist.events.keys() if type(key) is str])
        if gt_simulations != '':
            valid_stanames = sorted([key for key in self.statlist.stations.keys() if type(key) is str])
            valid_gt_components = ['CX','CY','CZ']

        # handle string displ_simulations
        if displ_simulations != '':
            if displ_simulations == 'all':
                self.all_tasks += [('displ',evid) for evid in valid_evids]
            elif displ_simulations.startswith('all-except:'):
                given_evids = displ_simulations[11:].split(',')
                given_invalid_evids = [evid for evid in given_evids if not evid in valid_evids]
                if len(given_invalid_evids) > 0:
                    self.log("### STOP : the following invalid eventIDs (that are not contained in FILE_EVENT_LIST) "+
                             "were detected in string 'displ_simulations':\n"+
                             "### "+',  '.join(given_invalid_evids)+"\n\n")
                    raise Exception("there were "+int(len(given_invalid_evids))+
                                    " invalid eventIDs detected in string 'displ_simulations'; see logfile '"+logfile+"'")
                self.all_tasks += [('displ',evid) for evid in valid_evids if not evid in given_evids]
            else:
                # expect here simply a ','-separated list of eventIDs
                given_evids = displ_simulations.split(',')
                given_invalid_evids = [evid for evid in given_evids if not evid in valid_evids]
                if len(given_invalid_evids) > 0:
                    self.log("### STOP : the following invalid eventIDs (that are not contained in FILE_EVENT_LIST) "+
                             "were detected in string 'displ_simulations':\n"+
                             "### "+',  '.join(given_invalid_evids)+"\n\n")
                    raise Exception("there were "+int(len(given_invalid_evids))+
                                    " invalid eventIDs detected in string 'displ_simulations'; see logfile '"+logfile+"'")
                self.all_tasks += [('displ',evid) for evid in given_evids]

        # handle string gt_simulations
        if gt_simulations != '':
            # first check gt_components in case of gt_simulations being "all" or starting with "all-except:"
            if gt_simulations == 'all' or gt_simulations.startswith('all-except:'):
                given_components = gt_components.split(',')
                if len(given_components) == 0:
                    self.log("### STOP : the string 'gt_components' is empty! In case of 'gt_simulations = all' "+
                             "or gt_simulations starting with 'all-except:', "+
                             "'gt_components' must contain at least one valid component\n\n")
                    raise Exception("no gt_components present in case of 'gt_simulations = all' or gt_simulations starting with 'all-except:' ; see logfile '"+logfile+"'")
                given_invalid_components = [comp for comp in given_components if not comp in valid_gt_components]
                if len(given_invalid_components) > 0:
                    self.log("### STOP : the following invalid components were detected in string 'gt_components':\n"+
                             "### "+',  '.join(given_invalid_components)+"\n"+
                             "### currently supported components are:\n### "+', ',join(valid_gt_components)+"\n\n")
                    raise Exception("there were "+int(len(given_invalid_components))+
                                    " invalid components detected in string 'gt_components'; see logfile '"+logfile+"'")
            # now handle 'gt_simulations = all'
            if gt_simulations == 'all':
                self.all_tasks += [('gt',staname+'_'+comp) for staname in valid_stanames 
                                   for comp in given_components]
                file_content = str(len(given_components))+"\n"+"\n".join(given_components)+"\n"
                self.gt_comp_files += [(staname,file_content) for staname in valid_stanames]
            # now handle 'gt_simulations = all-except:...'
            elif gt_simulations.startswith('all-except:'):
                given_stanames = gt_simulations[11:].split(',')
                given_invalid_stanames = [staname for staname in given_stanames if not staname in valid_stanames]
                if len(given_invalid_stanames) > 0:
                    self.log("### STOP : the following invalid station names (that are not contained in FILE_STATION_LIST) "+
                             "were detected in string 'gt_simulations':\n"+
                             "### "+',  '.join(given_invalid_stanames)+"\n\n")
                    raise Exception("there were "+int(len(given_invalid_stanames))+
                                    " invalid station names detected in string 'gt_simulations'; see logfile '"+logfile+"'")
                self.all_tasks += [('gt',staname+'_'+comp) for staname in valid_stanames 
                                   if not staname in given_stanames
                                   for comp in given_components]
                file_content = str(len(given_components))+"\n"+"\n".join(given_components)+"\n"
                self.gt_comp_files += [(staname,file_content) for staname in valid_stanames 
                                       if not staname in given_stanames]
            # now handle 'gt_simulations = specific'
            elif gt_simulations == 'specific':
                given_stanames_comps = gt_components.split(';')
                if len(given_stanames_comps) == 0:
                    self.log("### STOP : the string 'gt_components' is empty! In case of 'gt_simulations = specific', "+
                             "'gt_components' must contain a list of valid station name - component combinations\n\n")
                    raise Exception("no gt_components specification present in case of 'gt_simulations = specific' ; see logfile '"+logfile+"'")
                for i,staname_comps in enumerate(given_stanames_comps):
                    staname_comps_split = staname_comps.split(':')
                    if len(staname_comps_split) != 2:
                        self.log("### STOP : the "+str(i+1)+"'th station-component definition '"+staname_comps+
                                 "' of string 'gt_components' is invalid! 'gt_components' must be of form "+
                                 " 'Geophone23:CX,CY,CZ;Geophone2:CZ;ReceiverB:CY,CX' in case of 'gt_simulations = specific'\n\n")
                        raise Exception("gt_components specification invalid in case of 'gt_simulations = specific' ; see logfile '"+logfile+"'")
                    staname = staname_comps_split[0]
                    if not staname in valid_stanames:
                        self.log("### STOP : the "+str(i+1)+"'th station '"+staname+
                                 "' of string 'gt_components' is not in stations list! 'gt_components' must be of form "+
                                 " 'Geophone23:CX,CY,CZ;Geophone2:CZ;ReceiverB:CY,CX' in case of 'gt_simulations = specific'\n\n")
                        raise Exception("gt_components specification invalid in case of 'gt_simulations = specific' ; see logfile '"+logfile+"'")
                    given_components = staname_comps_split[1].split(',')
                    if len(given_components) == 0:
                        self.log("### STOP : the "+str(i+1)+"'th list of components '"+staname_comps_split[1]+
                                 "' of string 'gt_components' does not contain any components! 'gt_components' must be of form "+
                                 " 'Geophone23:CX,CY,CZ;Geophone2:CZ;ReceiverB:CY,CX' in case of 'gt_simulations = specific'\n\n")
                        raise Exception("gt_components specification invalid in case of 'gt_simulations = specific' ; see logfile '"+logfile+"'")
                    given_invalid_components = [comp for comp in given_components if not comp in valid_gt_components]
                    if len(given_invalid_components) > 0:
                        self.log("### STOP : the following invalid components were detected in the "+str(i+1)+
                                 "'th list of components '"+staname_comps_split[1]+"' of string 'gt_components':\n"+
                                 "### "+',  '.join(given_invalid_components)+"\n"+
                                 "### currently supported components are:\n"+
                                 "### "+', ',join(valid_gt_components)+
                                 "### 'gt_components' must be of form 'Geophone23:CX,CY,CZ;Geophone2:CZ;ReceiverB:CY,CX' in case of 'gt_simulations = specific'\n\n")
                        raise Exception("gt_components specification invalid in case of 'gt_simulations = specific' ; see logfile '"+logfile+"'")
                    # after all checks, it was verified that staname is a valid station name and that 
                    # all given_components for this station are valid. so add those Green functions to tasks list
                    self.all_tasks += [('gt',staname+'_'+comp) for comp in given_components]
                    file_content = str(len(given_components))+"\n"+"\n".join(given_components)+"\n"
                    self.gt_comp_files += [(staname,file_content)]
            # if gt_silmulations is neiter 'all', 'specific', nor starts with 'all-except:', raise an ERROR
            else:
                self.log("### STOP : gt_simulations has the value '"+gt_simulations+"'. It can be either equal to 'all' or "+
                         "'specific' or can start with 'all-except:'\n\n")
                raise Exception("invalid value of string 'gt_simulations'; see logfile '"+logfile+"'")

        # handle string measured_data_simulations
        if measured_data_simulations != '':
            if measured_data_simulations == 'all':
                self.all_tasks += [('data',evid) for evid in valid_evids]
            elif measured_data_simulations.startswith('all-except:'):
                given_evids = measured_data_simulations[11:].split(',')
                given_invalid_evids = [evid for evid in given_evids if not evid in valid_evids]
                if len(given_invalid_evids) > 0:
                    self.log("### STOP : the following invalid eventIDs (that are not contained in FILE_EVENT_LIST) "+
                             "were detected in string 'measured_data_simulations':\n"+
                             "### "+',  '.join(given_invalid_evids)+"\n\n")
                    raise Exception("there were "+int(len(given_invalid_evids))+
                                    " invalid eventIDs detected in string 'measured_data_simulations'; see logfile '"+logfile+"'")
                self.all_tasks += [('data',evid) for evid in valid_evids if not evid in given_evids]
            else:
                # expect here simply a ','-separated list of eventIDs
                given_evids = measured_data_simulations.split(',')
                given_invalid_evids = [evid for evid in given_evids if not evid in valid_evids]
                if len(given_invalid_evids) > 0:
                    self.log("### STOP : the following invalid eventIDs (that are not contained in FILE_EVENT_LIST) "+
                             "were detected in string 'measured_data_simulations':\n"+
                             "### "+',  '.join(given_invalid_evids)+"\n\n")
                    raise Exception("there were "+int(len(given_invalid_evids))+
                                    " invalid eventIDs detected in string 'measured_data_simulations'; see logfile '"+logfile+"'")
                self.all_tasks += [('data',evid) for evid in given_evids]
#
#-----------------------------------------------------------
#
    def iterate(self):
        for i,typ_id in enumerate(self.all_tasks):
            #
            # get current time, start of this iteration
            t_begin_iteration = time_time()
            #
            typ,sid = typ_id
                
            self.log('################################################################################\n'+
                     time_ctime()+' -- now doing the '+str(i+1)+'-th simulation out of '+str(len(self.all_tasks))+'\n'+
                     "TYPE '"+typ+"', ID '"+sid+"'\n"+
                     "\n")
            #
            # set all relevant parameter files correctly for this simulation
            self.log("set SPECFEM3D_Cartesian parameter files for this simulation now\n")
            try:
                self.setSpecfemCartParameters(typ,sid)
            except:
                self.log("### STOP : an error occurred setting the SPECFEM3D_Cartesian parameter files for this simulation\n")
                raise
            self.log("done\n"+
                     "\n")
            #
            # call some command by system call which conducts the simulation
            #
            # according to flag use_different_command_in_first_simulation, select which run command to use here
            if use_different_command_in_first_simulation and i==0:
                run_command = command_system_call_first_simulation
            else:
                run_command = command_system_call
            self.log("run command '"+run_command+"' now via system call\n")
            os_system(run_command)
            self.check_if_simulation_was_successful(typ)
            self.log("done\n\n")
            #
            # move OUTPUT_FILES
            self.log("copy SPECFEM3D_Cartesian OUTPUT_FILES\n")
            try:
                self.copySpecfemCartOutput()
            except:
                self.log("### STOP : an error occurred moving the standard SPECFEM3D_Cartesian output files\n")
                raise
            self.log("done\n\n")
            #
            # get current time, end of this iteration
            t_end_iteration = time_time()
            #
            # elapsed time for this iteration in h,min,sec
            this_min,this_sec = divmod(int(t_end_iteration-t_begin_iteration),60)
            this_h,this_min = divmod(this_min,60)
            #
            # mean elapsed time per simulation (computed from all iterations done, so far)
            t_mean_sec = (t_end_iteration-self.time_start)/(i+1)
            mean_min,mean_sec = divmod(int(t_mean_sec),60)
            mean_h,mean_min = divmod(mean_min,60)
            #
            # presumable time to finish
            t_finish = self.time_start+len(self.all_tasks)*t_mean_sec
            #
            # log all this time statistics
            self.log("current time -- "+time_ctime(t_end_iteration)+"\n"+
                     "elapsed time for this simulation (h:min:sec) -- %i : %i : %i"%(this_h,this_min,this_sec)+"\n"
                     "current mean elapsed time per simulation after %i simulations (h:min:sec) -- %i : %i : %i"%(i+1,mean_h,mean_min,mean_sec)+"\n"
                     "assuming all remaining "+str(len(self.all_tasks)-i-1)+" simulations to last the same time, script will finish presumably -- "+time_ctime(t_finish)+"\n\n\n")
            #
            # send logfile via email
            if send_emails and i in self.index_iteration_send_email:
                self.email_log("done %i out of %i; finish " % (i+1,len(self.all_tasks)) + time_ctime(t_finish))
            #
            # copy logfile to current SPECFEM3D_Cartesian OUTPUT_FILES directory
            # which was created in call self.copySpecfemCartOutput() above
            os_system('cp '+logfile+' '+self.outfile_base+'_OUTPUT_FILES')
#
#-----------------------------------------------------------
#
    def log(self,message):
        open(logfile,'a').write(message)
#
#-----------------------------------------------------------
#
    def email_log(self,subject):
        if email_receiver != '':
            os_system('mail -s "'+subject+'" '+email_receiver+' < '+logfile+' -- -f '+email_sender)
#
#-----------------------------------------------------------
#
    def setSpecfemCartParameters(self,typ,sid):
        # read in the parameter file DATA/Par_file_ASKI
        try:
            Par_file_ASKI = inputParameter(os_path.join(IN_DATA_FILES_PATH,'Par_file_ASKI'))
        except:
            self.log("   ERROR! could not create inputParameter object for Par_file_ASKI '"+
                     os_path.join(IN_DATA_FILES_PATH,'Par_file_ASKI')+"\n")
            raise 
        noKeys = Par_file_ASKI.keysNotPresent(['OVERWRITE_ASKI_OUTPUT'])
        if len(noKeys) > 0:
            self.log("   ERROR! in setSpecfemCartParameters: the following keywords are required in Par_file_ASKI '"+
                     os_path.join(IN_DATA_FILES_PATH,'Par_file_ASKI')+"':\n"+
                     "   "+',  '.join(noKeys)+"\n")
            raise Exception("missing keywords in Par_file_ASKI; see logfile '"+logfile+"'")
        #
        # remember current value of 'OVERWRITE_ASKI_OUTPUT' for this simulation
        self.overwrite_ASKI_output = Par_file_ASKI.lval('OVERWRITE_ASKI_OUTPUT')
        #
        # read in the parameter file DATA/Par_file
        try:
            Par_file = inputParameter(os_path.join(IN_DATA_FILES_PATH,'Par_file'))
        except:
            self.log("   ERROR! could not create inputParameter object for Par_file '"+
                     os_path.join(IN_DATA_FILES_PATH,'Par_file')+"\n")
            raise 
        noKeys = Par_file.keysNotPresent(['DT'])
        if len(noKeys) > 0:
            self.log("   ERROR! in setSpecfemCartParameters: the following keywords are required in Par_file '"+
                     os_path.join(IN_DATA_FILES_PATH,'Par_file')+"':\n"+
                     "   "+',  '.join(noKeys)+"\n")
            raise Exception("missing keywords in Par_file; see logfile '"+logfile+"'")
        #
        # remember current value of 'DT' for this simulation, accounting for any possible fortran notation for exponentials
        self.DT = float(Par_file.sval('DT').replace('D','e').replace('d','e'))
        #
        # now, according to simulation type, set parameter files
        #
        ##########################
        # typ=='displ':
        ##########################
        if typ=='displ':
            slat = self.evlist.events[sid]['slat']
            slon = self.evlist.events[sid]['slon']
            sdepth = self.evlist.events[sid]['sdepth']
            styp = self.evlist.events[sid]['styp']
            df = self.mparam.sval('MEASURED_DATA_FREQUENCY_STEP')
            nf = self.iparam.sval('ITERATION_STEP_NUMBER_OF_FREQ')
            jf = self.iparam.sval('ITERATION_STEP_INDEX_OF_FREQ')
            self.outfile_base = os_path.join(self.iter_path,self.iparam.sval('PATH_KERNEL_DISPLACEMENTS'),'kernel_displ_'+sid)
            self.log("   in setSpecfemCartParameters:\n"+
                     "      this is event evid = "+sid+"\n"+
                     "      source lat (X), lon (Y), depth (Z) = "+', '.join([slat,slon,sdepth])+"\n"+
                     "      source type = "+styp+"\n")
            if int(styp) == 0 and self.evlist.events[sid].has_key('force'):
                smag = self.evlist.events[sid]['mag']
                force = self.evlist.events[sid]['force']
                self.log("      force vector =  "+'  '.join(force)+"\n"+
                         "      force factor =  "+smag+"\n"+
                         "      nf,df = "+', '.join([nf,df])+"\n"+
                         "      frequency indices = "+jf+"\n"+
                         "      kernel displacement output file (basename) = '"+self.outfile_base+"'\n")
                Fx,Fy,Fz = force
                USE_FORCE_POINT_SOURCE = '.true.'
                try:
                    # BE AWARE! THAT IN SPECFEM THE FIRST (wavefield point) COORDINATE (x) IS LON, AND THE SECOND COORDINATE (y) IS LAT
                    # IN CARTESIAN ASKI APPLICATIONS, HOWEVER, THE FIRST (wavefield point) COORDINATE IS THE FIRST COLUMN ( = LAT) IN 
                    # EVENT- AND STATION LIST FILES, AND THE SECOND COORDINATE (column) IS LON, SO INTERCHANGE LAT AND LON HERE!!!
                    setForcesolution(hdur=str(5*self.DT),lat=slon,lon=slat,depth=sdepth,factor_force=smag,FE=Fx,FN=Fy,FUP=Fz)
                except:
                    self.log("   ERROR! in setSpecfemCartParameters: exception raised while setting FORCESOLUTION\n")
                    raise
            elif int(styp) == 1 and self.evlist.events[sid].has_key('momten'):
                momten = self.evlist.events[sid]['momten']
                momten_DynCm = Moment_tensor_Nm2DynCm(momten)
                self.log("      moment tensor in Nm =  "+'  '.join(momten)+" (read from ASKI event list file)\n"+
                         "      moment tensor in dyn*cm =  "+'  '.join(momten_DynCm)+" (write to SPECFEM CMTSOLUTION file)\n"+
                         "      nf,df = "+', '.join([nf,df])+"\n"+
                         "      frequency indices = "+jf+"\n"+
                         "      kernel displacement output file (basename) = '"+self.outfile_base+"'\n")
                Mrr,Mtt,Mpp,Mrt,Mrp,Mtp = momten_DynCm
                USE_FORCE_POINT_SOURCE = '.false.'
                try:
                    # BE AWARE! THAT IN SPECFEM THE FIRST (wavefield point) COORDINATE (x) IS LON, AND THE SECOND COORDINATE (y) IS LAT
                    # IN CARTESIAN ASKI APPLICATIONS, HOWEVER, THE FIRST (wavefield point) COORDINATE IS THE FIRST COLUMN ( = LAT) IN 
                    # EVENT- AND STATION LIST FILES, AND THE SECOND COORDINATE (column) IS LON, SO INTERCHANGE LAT AND LON HERE!!!
                    setCmtsolution(evname=sid,hdur=str(5*self.DT),lat=slon,lon=slat,depth=sdepth,Mrr=Mrr,Mtt=Mtt,Mpp=Mpp,Mrt=Mrt,Mrp=Mrp,Mtp=Mtp)
                except:
                    self.log("   ERROR! in setSpecfemCartParameters: exception raised while setting CMTSOLUTION\n")
                    raise
            else:
                self.log("   ERROR! in setSpecfemCartParameters: source type is undefined (neither force nor moment tensor)\n")
                raise Exception("source type is undefined")

            if create_specfem_stations:
                log_string = ("      the SPECFEM STATIONS file was created from the "+str(self.statlist.nstat)+" stations in ASKI's station list file")
                # BE AWARE! THAT IN SPECFEM THE FIRST (wavefield point) COORDINATE (x) IS LON, AND THE SECOND COORDINATE (y) IS LAT
                # IN CARTESIAN ASKI APPLICATIONS, HOWEVER, THE FIRST (wavefield point) COORDINATE IS THE FIRST COLUMN ( = LAT) IN 
                # EVENT- AND STATION LIST FILES, AND THE SECOND COORDINATE (column) IS LON, SO INTERCHANGE LAT AND LON HERE!!!
                STATIONS_content = '\n'.join(['   '.join([self.statlist.stations[i]['staname'], self.statlist.stations[i]['netcode'], 
                                                          self.statlist.stations[i]['lon'], self.statlist.stations[i]['lat'],
                                                          '0.0', self.statlist.stations[i]['alt']
                                                          ])
                                              for i in range(self.statlist.nstat)
                                              ])+'\n'
                try:
                    open(os_path.join(IN_DATA_FILES_PATH,'STATIONS'),'w').write(STATIONS_content)
                except:
                    self.log("   ERROR! in setSpecfemCartParameters: could not open STATIONS file '"+os_path.join(IN_DATA_FILES_PATH,'STATIONS')+
                                    "' to write\n")
                    raise
                self.log(log_string+"\n")
            else:
                self.log("      the SPECFEM STATIONS file was not modified by this script, you should have set it yourself correctly\n")

            if os_path.exists(self.outfile_base+'_OUTPUT_FILES'):
                if self.overwrite_ASKI_output:
                    self.log("      SPECFEM3D_Cartesian OUTPUT_FILES will be copied to '"+
                             self.outfile_base+'_OUTPUT_FILES'+"' (exists and will be overwritten).\n")
                else:
                    self.log("   ERROR! in setSpecfemCartParameters: path for SPECFEM3D_Cartesian OUTPUT_FILES '"+
                             self.outfile_base+'_OUTPUT_FILES'+"' exists and must not be overwritten "+
                             "(according to 'OVERWRITE_ASKI_OUTPUT' in Par_file_ASKI)\n")
                    raise Exception('path for SPECFEM3D_Cartesian OUTPUT_FILES already exists')
            else:
                self.log("      SPECFEM3D_Cartesian OUTPUT_FILES will be copied to '"+
                         self.outfile_base+'_OUTPUT_FILES'+"'.\n")

            try:
                setParfile(os_path.join(IN_DATA_FILES_PATH,'Par_file'),
                           [('USE_FORCE_POINT_SOURCE',USE_FORCE_POINT_SOURCE),
                            ('USE_RICKER_TIME_FUNCTION','.false.')])
            except:
                self.log("   ERROR! in setSpecfemCartParameters: exception raised while setting Par_file\n")
                raise
            try:
                setParfile(os_path.join(IN_DATA_FILES_PATH,'Par_file_ASKI'),
                           [('COMPUTE_ASKI_OUTPUT','.true.'),('ASKI_outfile',self.outfile_base),
                            ('ASKI_output_ID',sid),('ASKI_df',df),('ASKI_nf',nf),('ASKI_jf',jf)])
                if define_ASKI_output_volume_by_inversion_grid:
                    setParfile(os_path.join(IN_DATA_FILES_PATH,'Par_file_ASKI'),
                               [('ASKI_type_inversion_grid',self.ASKI_type_inversion_grid),('ASKI_wx',self.ASKI_wx),
                                ('ASKI_wy',self.ASKI_wy),('ASKI_wz',self.ASKI_wz),('ASKI_rot_X',self.ASKI_rot_X),
                                ('ASKI_rot_Y',self.ASKI_rot_Y),('ASKI_rot_Z',self.ASKI_rot_Z),('ASKI_cx',self.ASKI_cx),
                                ('ASKI_cy',self.ASKI_cy),('ASKI_cz',self.ASKI_cz)])
            except:
                self.log("   ERROR! in setSpecfemCartParameters: exception raised while setting Par_file_ASKI\n")
                raise

        ##########################
        # typ=='gt':
        ##########################
        elif typ=='gt':
            staname_comp = sid.split('_')
            if len(staname_comp) != 2:
                self.log("   ERROR! in setSpecfemCartParameters: simulation ID '"+sid+"' is invalid for simulation "+
                         "type 'gt'. Must be of form 'stationname_component', where stationname MUST NOT contain '_' "+
                         "characters!\n")
                raise Exception("invalid simulation ID '"+sid+"' for Green tensor simulation")
            staname = staname_comp[0]
            comp = staname_comp[1]
            Fx=Fy=Fz='0.d0'
            if comp == 'CX':
                Fx='1.d0'
            elif comp == 'CY':
                Fy='1.d0'
            elif comp == 'CZ':
                Fz='1.d0'
            nwname = self.statlist.stations[staname]['netcode']
            lat = self.statlist.stations[staname]['lat']
            lon = self.statlist.stations[staname]['lon']
            # use value of 'alt' to put on CMTSOLUTION line 'depth:', in order to allow receivers not to be located on the surface
            # it makes sense to use 'USE_SOURCES_RECVS_Z = .true.' in constants.h, as then depth and , hence, the value for 'alt'
            # is interpreted as the 'Z' coordinate
            alt = self.statlist.stations[staname]['alt']
            df = self.mparam.sval('MEASURED_DATA_FREQUENCY_STEP')
            nf = self.iparam.sval('ITERATION_STEP_NUMBER_OF_FREQ')
            jf = self.iparam.sval('ITERATION_STEP_INDEX_OF_FREQ')
            self.outfile_base = os_path.join(self.iter_path,self.iparam.sval('PATH_KERNEL_GREEN_TENSORS'),'kernel_gt_'+sid)
            self.log("   in setSpecfemCartParameters:\n"+
                     "      this is station stname,nwname = "+', '.join([staname,nwname])+"\n"+
                     "      lat (X), lon (Y), alt (Z, used as depth value in FORCESOLUTION) = "+', '.join([lat,lon,alt])+"\n"+
                     "      green tensor component = "+comp+"\n"+
                     "      nf,df = "+', '.join([nf,df])+"\n"+
                     "      frequency indices = "+jf+"\n"+
                     "      kernel green tensor output file (basename) = '"+self.outfile_base+"'\n")

            if create_specfem_stations:
                log_string = ("      the SPECFEM STATIONS file was created from the "+str(self.statlist.nstat)+" stations in ASKI's station list file")
                # BE AWARE! THAT IN SPECFEM THE FIRST (wavefield point) COORDINATE (x) IS LON, AND THE SECOND COORDINATE (y) IS LAT
                # IN CARTESIAN ASKI APPLICATIONS, HOWEVER, THE FIRST (wavefield point) COORDINATE IS THE FIRST COLUMN ( = LAT) IN 
                # EVENT- AND STATION LIST FILES, AND THE SECOND COORDINATE (column) IS LON, SO INTERCHANGE LAT AND LON HERE!!!
                STATIONS_content = '\n'.join(['   '.join([self.statlist.stations[i]['staname'], self.statlist.stations[i]['netcode'], 
                                                          self.statlist.stations[i]['lon'], self.statlist.stations[i]['lat'],
                                                          '0.0', self.statlist.stations[i]['alt']
                                                          ])
                                              for i in range(self.statlist.nstat)
                                              ])+'\n'
                try:
                    open(os_path.join(IN_DATA_FILES_PATH,'STATIONS'),'w').write(STATIONS_content)
                except:
                    self.log("   ERROR! in setSpecfemCartParameters: could not open STATIONS file '"+os_path.join(IN_DATA_FILES_PATH,'STATIONS')+
                                    "' to write\n")
                    raise
                self.log(log_string+"\n")
            else:
                self.log("      the SPECFEM STATIONS file was not modified by this script, you should have set it yourself correctly\n")

            if os_path.exists(self.outfile_base+'_OUTPUT_FILES'):
                if self.overwrite_ASKI_output:
                    self.log("      SPECFEM3D_Cartesian OUTPUT_FILES will be copied to '"+
                             self.outfile_base+'_OUTPUT_FILES'+"' (exists and will be overwritten).\n")
                else:
                    self.log("   ERROR! in setSpecfemCartParameters: path for SPECFEM3D_Cartesian OUTPUT_FILES '"+
                             self.outfile_base+'_OUTPUT_FILES'+"' exists and must not be overwritten "+
                             "(according to 'OVERWRITE_ASKI_OUTPUT' in Par_file_ASKI)\n")
                    raise Exception('path for SPECFEM3D_Cartesian OUTPUT_FILES already exists')
            else:
                self.log("      SPECFEM3D_Cartesian OUTPUT_FILES will be copied to '"+
                         self.outfile_base+'_OUTPUT_FILES'+"'.\n")
            try:
                # BE AWARE! THAT IN SPECFEM THE FIRST (wavefield point) COORDINATE (x) IS LON, AND THE SECOND COORDINATE (y) IS LAT
                # IN CARTESIAN ASKI APPLICATIONS, HOWEVER, THE FIRST (wavefield point) COORDINATE IS THE FIRST COLUMN ( = LAT) IN 
                # EVENT- AND STATION LIST FILES, AND THE SECOND COORDINATE (column) IS LON, SO INTERCHANGE LAT AND LON HERE!!!
                setForcesolution(hdur=str(5*self.DT),lat=lon,lon=lat,depth=alt,factor_force='1.d0',FE=Fx,FN=Fy,FUP=Fz)
            except:
                self.log("   ERROR! in setSpecfemCartParameters: exception raised while setting FORCESOLUTION\n")
                raise
            try:
                setParfile(os_path.join(IN_DATA_FILES_PATH,'Par_file'),
                           [('USE_FORCE_POINT_SOURCE','.true.'),('USE_RICKER_TIME_FUNCTION','.false.')])
            except:
                self.log("   ERROR! in setSpecfemCartParameters: exception raised while setting Par_file\n")
                raise
            try:
                setParfile(os_path.join(IN_DATA_FILES_PATH,'Par_file_ASKI'),
                           [('COMPUTE_ASKI_OUTPUT','.true.'),('ASKI_outfile',self.outfile_base),
                            ('ASKI_output_ID',sid),('ASKI_df',df),('ASKI_nf',nf),('ASKI_jf',jf)])
                if define_ASKI_output_volume_by_inversion_grid:
                    setParfile(os_path.join(IN_DATA_FILES_PATH,'Par_file_ASKI'),
                               [('ASKI_type_inversion_grid',self.ASKI_type_inversion_grid),('ASKI_wx',self.ASKI_wx),
                                ('ASKI_wy',self.ASKI_wy),('ASKI_wz',self.ASKI_wz),('ASKI_rot_X',self.ASKI_rot_X),
                                ('ASKI_rot_Y',self.ASKI_rot_Y),('ASKI_rot_Z',self.ASKI_rot_Z),('ASKI_cx',self.ASKI_cx),
                                ('ASKI_cy',self.ASKI_cy),('ASKI_cz',self.ASKI_cz)])
            except:
                self.log("   ERROR! in setSpecfemCartParameters: exception raised while setting Par_file_ASKI\n")
                raise

        ##########################
        # typ=='data':
        ##########################
        elif typ=='data':
            slat = self.evlist.events[sid]['slat']
            slon = self.evlist.events[sid]['slon']
            sdepth = self.evlist.events[sid]['sdepth']
            styp = self.evlist.events[sid]['styp']
            self.outfile_base = (self.mparam.sval('PATH_MEASURED_DATA')+'data_'+sid)
            self.log("   in setSpecfemCartParameters:\n"+
                     "      this is event evid = "+sid+"\n"+
                     "      source source lat (X), lon (Y), depth (Z) = "+', '.join([slat,slon,sdepth])+"\n"+
                     "      source type = "+styp+"\n")
            if int(styp) == 0 and self.evlist.events[sid].has_key('force'):
                smag = self.evlist.events[sid]['mag']
                force = self.evlist.events[sid]['force']
                self.log("      force vector =  "+'  '.join(force)+"\n"+
                         "      force factor =  "+smag+"\n")
                Fx,Fy,Fz = force
                USE_FORCE_POINT_SOURCE = '.true.'
                try:
                    # in case of data-simulations, do not touch hdur (so, if set manually previously, 
                    # you can use any source time function)
                    # BE AWARE, THAT IN CASE OF MIXED data/displ/gt SIMULATIONS, hdur WILL HAVE CHANGED
                    # IN gt/displ SIMULATIONS TO 5*self.DT !!!
                    #
                    # BE AWARE! THAT IN SPECFEM THE FIRST (wavefield point) COORDINATE (x) IS LON, AND THE SECOND COORDINATE (y) IS LAT
                    # IN CARTESIAN ASKI APPLICATIONS, HOWEVER, THE FIRST (wavefield point) COORDINATE IS THE FIRST COLUMN ( = LAT) IN 
                    # EVENT- AND STATION LIST FILES, AND THE SECOND COORDINATE (column) IS LON, SO INTERCHANGE LAT AND LON HERE!!!
                    setForcesolution(lat=slon,lon=slat,depth=sdepth,factor_force=smag,FE=Fx,FN=Fy,FUP=Fz)
                except:
                    self.log("   ERROR! in setSpecfemCartParameters: exception raised while setting FORCESOLUTION\n")
                    raise
            elif int(styp) == 1 and self.evlist.events[sid].has_key('momten'):
                momten = self.evlist.events[sid]['momten']
                self.log("      moment tensor in Nm =  "+'  '.join(momten)+" (read from ASKI event list file)\n")
                momten_DynCm = Moment_tensor_Nm2DynCm(momten)
                self.log("      moment tensor in dyn*cm =  "+'  '.join(momten_DynCm)+" (write to SPECFEM CMTSOLUTION file)\n")
                Mrr,Mtt,Mpp,Mrt,Mrp,Mtp = momten_DynCm
                try:
                    # BE AWARE! THAT IN SPECFEM THE FIRST (wavefield point) COORDINATE (x) IS LON, AND THE SECOND COORDINATE (y) IS LAT
                    # IN CARTESIAN ASKI APPLICATIONS, HOWEVER, THE FIRST (wavefield point) COORDINATE IS THE FIRST COLUMN ( = LAT) IN 
                    # EVENT- AND STATION LIST FILES, AND THE SECOND COORDINATE (column) IS LON, SO INTERCHANGE LAT AND LON HERE!!!
                    setCmtsolution(evname=sid,lat=slon,lon=slat,depth=sdepth,Mrr=Mrr,Mtt=Mtt,Mpp=Mpp,Mrt=Mrt,Mrp=Mrp,Mtp=Mtp)
                except:
                    self.log("   ERROR! in setSpecfemCartParameters: exception raised while setting CMTSOLUTION\n")
                    raise
            else:
                self.log("   ERROR! in setSpecfemCartParameters: source type is undefined (neither force nor moment tensor)\n")
                raise Exception("source type is undefined")

            if create_specfem_stations:
                log_string = ("      the SPECFEM STATIONS file was created from the "+str(self.statlist.nstat)+" stations in ASKI's station list file")
                # BE AWARE! THAT IN SPECFEM THE FIRST (wavefield point) COORDINATE (x) IS LON, AND THE SECOND COORDINATE (y) IS LAT
                # IN CARTESIAN ASKI APPLICATIONS, HOWEVER, THE FIRST (wavefield point) COORDINATE IS THE FIRST COLUMN ( = LAT) IN 
                # EVENT- AND STATION LIST FILES, AND THE SECOND COORDINATE (column) IS LON, SO INTERCHANGE LAT AND LON HERE!!!
                STATIONS_content = '\n'.join(['   '.join([self.statlist.stations[i]['staname'], self.statlist.stations[i]['netcode'], 
                                                          self.statlist.stations[i]['lon'], self.statlist.stations[i]['lat'],
                                                          '0.0', self.statlist.stations[i]['alt']
                                                          ])
                                              for i in range(self.statlist.nstat)
                                              ])+'\n'
                try:
                    open(os_path.join(IN_DATA_FILES_PATH,'STATIONS'),'w').write(STATIONS_content)
                except:
                    self.log("   ERROR! in setSpecfemCartParameters: could not open STATIONS file '"+os_path.join(IN_DATA_FILES_PATH,'STATIONS')+
                                    "' to write\n")
                    raise
                self.log(log_string+"\n")
            else:
                self.log("      the SPECFEM STATIONS file was not modified by this script, you should have set it yourself correctly\n")

            if os_path.exists(self.outfile_base+'_OUTPUT_FILES'):
                if self.overwrite_ASKI_output:
                    self.log("      SPECFEM3D_Cartesian OUTPUT_FILES will be copied to '"+
                             self.outfile_base+'_OUTPUT_FILES'+"' (exists and will be overwritten).\n")
                else:
                    self.log("   ERROR! in setSpecfemCartParameters: path for SPECFEM3D_Cartesian OUTPUT_FILES '"+
                             self.outfile_base+'_OUTPUT_FILES'+"' exists and must not be overwritten "+
                             "(according to 'OVERWRITE_ASKI_OUTPUT' in Par_file_ASKI)\n")
                    raise Exception('path for SPECFEM3D_Cartesian OUTPUT_FILES already exists')
            else:
                self.log("      SPECFEM3D_Cartesian OUTPUT_FILES will be copied to '"+
                         self.outfile_base+'_OUTPUT_FILES'+"'.\n")

            try:
                # in case of data-simulations, do not touch USE_RICKER_TIME_FUNCTION 
                # (so, if set manually previously, you can use any source time function)
                setParfile(os_path.join(IN_DATA_FILES_PATH,'Par_file'),[('USE_FORCE_POINT_SOURCE',USE_FORCE_POINT_SOURCE)])
            except:
                self.log("   ERROR! in setSpecfemCartParameters: exception raised while setting Par_file\n")
                raise
            try:
                setParfile(os_path.join(IN_DATA_FILES_PATH,'Par_file_ASKI'),[('COMPUTE_ASKI_OUTPUT','.false.')])
            except:
                self.log("   ERROR! in setSpecfemCartParameters: exception raised while setting Par_file_ASKI\n")
                raise
        else:
            self.log("   ERROR! in setSpecfemCartParameters: simulation type '"+typ+"' not supported\n")
            raise Exception("simulation type '"+typ+"' not supported")
#
#-----------------------------------------------------------
#
    def check_if_simulation_was_successful(self,typ):
        # check if all binaries exist which should have been compiled correctly
        bin_list = ['bin/xdecompose_mesh','bin/xgenerate_databases','bin/xspecfem3D']
        bin_not_exist = [f for f in bin_list if not os_path.exists(f)]
        if len(bin_not_exist)>0:
            self.log("### ERROR : the following binaries are not present, hence were not compiled correctly before starting this run:\n"+
                     "            '"+"', '".join(bin_not_exist)+"'\n")
            raise Exception("some binaries are not present; see logfile '"+logfile+"'")
        #
        # check if xdecompose_mesh_SCOTCH wrote all database files
        database_list = [os_path.join(LOCAL_PATH,'proc%6.6i_Database'%i) for i in range(self.nproc)]
        database_not_exist = [f for f in database_list if not os_path.exists(f)]
        if len(database_not_exist)>0:
            self.log("### ERROR : the following database files were not created by xdecompose_mesh: '"+
                     "', ".join(database_not_exist)+"'\n")
            raise Exception("some database files were not created by xdecompose_mesh; see logfile '"+
                            logfile+"'")
        #
        # check if database files are not empty
        database_empty = [f for f in database_list if os_path.getsize(f) <= 0]
        if len(database_empty)>0:
            self.log("### ERROR : the following database files are empty: '"+"', ".join(database_empty)+"'\n")
            raise Exception("some database files are empty; see logfile '"+logfile+"'")
        #
        # check if there are any "error_message*.txt" files in OUTPUT_FILES, i.e. if exit_mpi.f90 was called
        if any( [f.startswith('error_message') for f in os_listdir(OUTPUT_FILES_PATH)] ):
            self.log("### ERROR : there are files 'error_message*' in '"+OUTPUT_FILES_PATH+"'; "+
                     "mpiruns may have exited unintendedly\n")
            raise Exception("there are files 'error_message*' in '"+OUTPUT_FILES_PATH+"'; "+
                            "mpiruns may have exited unintendedly")
        #
        # check if there is a file "output_solver.txt"
        if not os_path.exists(os_path.join(OUTPUT_FILES_PATH,'output_solver.txt')):
            self.log("### ERROR : as there is no file 'output_solver.txt' in '"+OUTPUT_FILES_PATH+"', "+
                     "this suggests that the solver did not run correctly\n")
            raise Exception("as there is no file 'output_solver.txt' in '"+OUTPUT_FILES_PATH+"', "+
                     "this suggests that the solver did not run correctly")
        #
        # check if there is a file "starttimeloop.txt"
        if not os_path.exists(os_path.join(OUTPUT_FILES_PATH,'starttimeloop.txt')):
            self.log("### ERROR : as there is no file 'starttimeloop.txt' in '"+OUTPUT_FILES_PATH+"', "+
                     "this suggests that the solver did not run correctly\n")
            raise Exception("as there is no file 'starttimeloop.txt' in '"+OUTPUT_FILES_PATH+"', "+
                     "this suggests that the solver did not run correctly")
        #
        # check if one of the last lines of output_solver.txt is "End of the simulation"
        last_lines_output_solver = open(os_path.join(OUTPUT_FILES_PATH,'output_solver.txt'),'r').readlines()[-10:]
        if not any( ['End of the simulation' in line for line in last_lines_output_solver] ):
            self.log("### ERROR : the last 10 lines of file 'output_solver.txt' in '"+OUTPUT_FILES_PATH+"', "+
                     "do not contain line 'End of the simulation'\n")
            raise Exception("the last 10 lines of file 'output_solver.txt' in '"+OUTPUT_FILES_PATH+"', "+
                     "do not contain line 'End of the simulation'")
        #
        # check if in case of simulation types displ and gt, actually kernel output was produced successfully
        if typ in ['displ','gt']:
            if not os_path.exists(os_path.join(OUTPUT_FILES_PATH,'LOG_ASKI_finish.txt')):
                self.log("### ERROR : as there is no file 'LOG_ASKI_finish.txt' in '"+OUTPUT_FILES_PATH+"', "+
                         "although this is type '"+typ+"', this suggests that ASKI output was not produced correctly\n")
                raise Exception("as there is no file 'LOG_ASKI_finish.txt' in '"+OUTPUT_FILES_PATH+"', "+
                                "this suggests that ASKI output was not produced correctly")
            lines_ASKI_finish = open(os_path.join(OUTPUT_FILES_PATH,'LOG_ASKI_finish.txt'),'r').readlines()
            if not any( ["successfully created ASKI output, as specified in 'LOG_ASKI_start.txt'" in line for line in lines_ASKI_finish] ):
                self.log("### ERROR : file 'LOG_ASKI_finish.txt' in '"+OUTPUT_FILES_PATH+"', "+
                         "does not contain line 'successfully created ASKI output'\n")
                raise Exception("file 'LOG_ASKI_finish.txt' in '"+OUTPUT_FILES_PATH+"', "+
                                "does not contain line 'successfully created ASKI output'")
#
#-----------------------------------------------------------
#
    def copySpecfemCartOutput(self):
        if os_path.exists(self.outfile_base+'_OUTPUT_FILES'):
            if self.overwrite_ASKI_output:
                # remove existing output first, before copying new one
                self.log("   in copySpecfemCartOutput: calling \"os_system('rm -rf "+self.outfile_base+'_OUTPUT_FILES/*'+"')\"\n")
                os_system('rm -rf '+self.outfile_base+'_OUTPUT_FILES/*')
            #else: already raised Exception in this case in routine setSpecfemCartParameters
        else: 
            self.log("   in copySpecfemCartOutput: calling \"os_mkdir('"+self.outfile_base+'_OUTPUT_FILES'+"')\"\n")
            os_mkdir(self.outfile_base+'_OUTPUT_FILES')
        #
        # finally copy output files 
        self.log("   in copySpecfemCartOutput: calling \"os_system('cp "+os_path.join(OUTPUT_FILES_PATH,'*')+" "+
                 self.outfile_base+'_OUTPUT_FILES'+"')\"\n")
        os_system('cp '+os_path.join(OUTPUT_FILES_PATH,'*')+' '+self.outfile_base+'_OUTPUT_FILES')
############################################################
# END OF CLASS simulation
############################################################
#
#
#-----------------------------------------------------------
#
def Moment_tensor_Nm2DynCm(momten):
    #return [  'e+'.join([ m.lower().split('e+')[0] , str(int(m.lower().split('e+')[1])+7) ])  for m in momten  ]
    return [str(float(m)*1e+7) for m in momten]
#
#-----------------------------------------------------------
#
def setParfile(filename,keys_vals):
    # locally define dictionary for better handling of incoming key,value pairs which are to be changed in parameter file
    values = dict(keys_vals)
    keys = values.keys()
    #
    # read in parameter file
    try:
        lines_orig = open(filename,'r').readlines()
    except:
        raise Exception("could not open parameter file '"+filename+"' to read")
    #
    # iterate over all lines and modify line if necessary
    lines_new = []
    for line in lines_orig:
        # ignore empty lines, comment lines and invalid lines (which do not contain any '=' in front of a comment)
        if line.strip() == '' or line.strip()[0:1]=='#' or not '=' in line.split('#')[0]:
            lines_new.append(line)
            continue
        #
        # if the key,value pair of this valid line is to be modified, do so
        key_line = line.split('=')[0].strip()
        val_line = line.split('=')[1].split('#')[0].strip()
        if key_line in keys:
            # first of all, remove newline character from end of line and append it to modified line
            line = line.strip('\n')
            # replace old value by new value, but keep all whitespace and commentary of this line as was before
            key_part = line.split('=')[0]
            val_part = line.split('=')[1].split('#')[0]
            # comment part ends on newline character
            if '#' in line:
                comment_part = '#'+line.split('=')[1].split('#')[1]+'\n'
            else:
                comment_part = '\n'
            if val_line == '':
                val_part = ' '+values[key_line]+' '
            else:
                val_part = val_part.replace(val_line,values[key_line])
            line = key_part+'='+val_part+comment_part
            # remove key of this line from keys list, in order to check (in the end) if all keys were found in the file
            keys.remove(key_line)
        #
        # add the line (if modified or not) to list of new lines
        lines_new.append(line)
    #
    # check if there are any keys which were not found on valid lines in the file
    if len(keys)>0:
        raise Exception("could not find the following parameters in parameter file '"+filename+"': "+
                        "'"+"', '".join(keys)+"'")
    #
    # if every key was found and the respective value was modified, write modified lines to file
    try:
        open(filename,'w').writelines(lines_new)
    except:
        raise Exception("could not open parameter file '"+filename+"' to write modified lines")
#
#-----------------------------------------------------------
#
def setForcesolution(hdur=None,lat=None,lon=None,depth=None,factor_force=None,FE=None,FN=None,FUP=None):
    def setForcesolution_replaceLine(line,key,newval):
        line_split = line.split(key)
        if line_split[1].strip() == '':
            return line_split[0]+key+'   '+newval+'\n'
        else:
            return line_split[0]+key+(line_split[1].replace(line_split[1].strip(),newval))
    #
    # read lines of FORCESOLUTION file
    try:
        lines = open(os_path.join(IN_DATA_FILES_PATH,'FORCESOLUTION'),'r').readlines()
    except:
        raise Exception("could not open FORCESOLUTION file '"+os_path.join(IN_DATA_FILES_PATH,'FORCESOLUTION')+
                        "' to read")
    # now modify lines
    if hdur is not None:
        lines[2] = setForcesolution_replaceLine(lines[2],'f0:',hdur)
    if lat is not None:
        lines[3] = setForcesolution_replaceLine(lines[3],'latorUTM:',lat)
    if lon is not None:
        lines[4] = setForcesolution_replaceLine(lines[4],'longorUTM:',lon)
    if depth is not None:
        lines[5] = setForcesolution_replaceLine(lines[5],'depth:',depth)
    if factor_force is not None:
        lines[6] = setForcesolution_replaceLine(lines[6],'factor force source:',factor_force)
    if FE is not None:
        lines[7] = setForcesolution_replaceLine(lines[7],'component dir vect source E:',FE)
    if FN is not None:
        lines[8] = setForcesolution_replaceLine(lines[8],'component dir vect source N:',FN)
    if FUP is not None:
        lines[9] = setForcesolution_replaceLine(lines[9],'component dir vect source Z_UP:',FUP)
    # write modified lines to new FORCESOLUTION file
    try:
        open(os_path.join(IN_DATA_FILES_PATH,'FORCESOLUTION'),'w').writelines(lines)
    except:
        raise Exception("could not open FORCESOLUTION file '"+os_path.join(IN_DATA_FILES_PATH,'FORCESOLUTION')+
                        "' to write modified lines")
#
#-----------------------------------------------------------
#
def setCmtsolution(evname=None,hdur=None,lat=None,lon=None,depth=None,Mrr=None,Mtt=None,Mpp=None,Mrt=None,Mrp=None,Mtp=None):
    def setCmtsolution_replaceLine(line,key,newval):
        line_split = line.split(key)
        if line_split[1].strip() == '':
            return line_split[0]+key+'   '+newval+'\n'
        else:
            return line_split[0]+key+(line_split[1].replace(line_split[1].strip(),newval))
    #
    # read lines of CMTSOLUTION file
    try:
        lines = open(os_path.join(IN_DATA_FILES_PATH,'CMTSOLUTION'),'r').readlines()
    except:
        raise Exception("could not open CMTSOLUTION file '"+os_path.join(IN_DATA_FILES_PATH,'CMTSOLUTION')+
                        "' to read")
    # now modify lines
    if evname is not None:
        lines[1] = setCmtsolution_replaceLine(lines[1],'event name:',evname)
    if hdur is not None:
        lines[3] = setCmtsolution_replaceLine(lines[3],'half duration:',hdur)
    if lat is not None:
        lines[4] = setCmtsolution_replaceLine(lines[4],'latorUTM:',lat)
    if lon is not None:
        lines[5] = setCmtsolution_replaceLine(lines[5],'longorUTM:',lon)
    if depth is not None:
        lines[6] = setCmtsolution_replaceLine(lines[6],'depth:',depth)
    if Mrr is not None:
        lines[7] = setCmtsolution_replaceLine(lines[7],'Mrr:',Mrr)
    if Mtt is not None:
        lines[8] = setCmtsolution_replaceLine(lines[8],'Mtt:',Mtt)
    if Mpp is not None:
        lines[9] = setCmtsolution_replaceLine(lines[9],'Mpp:',Mpp)
    if Mrt is not None:
        lines[10] = setCmtsolution_replaceLine(lines[10],'Mrt:',Mrt)
    if Mrp is not None:
        lines[11] = setCmtsolution_replaceLine(lines[11],'Mrp:',Mrp)
    if Mtp is not None:
        lines[12] = setCmtsolution_replaceLine(lines[12],'Mtp:',Mtp)
    # write modified lines to new CMTSOLUTION file
    try:
        open(os_path.join(IN_DATA_FILES_PATH,'CMTSOLUTION'),'w').writelines(lines)
    except:
        raise Exception("could not open CMTSOLUTION file '"+os_path.join(IN_DATA_FILES_PATH,'CMTSOLUTION')+
                        "' to write modified lines")
#
############################################################
# MAIN
############################################################
#
def main():
    if runs_on_SGE:
        full_log_name = os_path.join(SGE_o_workdir,logfile)
    else:
        full_log_name = logfile
    open(logfile,'w').write('this is log '+full_log_name+'\nstarting script -- '+time_ctime()+'\n\n')

    # initiate simulation
    try:
        sm = simulation()
    except:
        open(logfile,'a').write('### STOP --'+time_ctime()+
                                '-- : could not initiate specfem3dCartesianForASKI simulations\n\n')
        raise

    # iterate over individual specfem3dCartesianForASKI simulations
    try:
        sm.iterate()
    except:
        open(logfile,'a').write('\n\n### STOP --'+time_ctime()+
                                '-- : there was an error iterating over the individual specfem3dCartesianForASKI simulations\n\n')
        if send_emails:
            sm.email_log('ERROR in one simulation')
        raise

    open(logfile,'a').write("\n"+
                            "--"+time_ctime()+"-- exiting script now. Good Bye\n"+
                            "\n")
    if send_emails:
        sm.email_log("successfully finished")
    sys_exit()
#
############################################################
#
if __name__ == "__main__":
    main()
