#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.config Configuration settings and paths for the eagle package.
#
# This module defines configuration settings and paths for the eagle package on various computing systems.
# The settings for each system are bundled in a dictionary with an appropriate name. In addition, the settings
# corresponding to the current system are copied into variables in the module's outermost namespace.
# For example, after executing "import eagle.config" the path to the SKIRT executable for
# the system on which the script is running can be found as the value of "eagle.config.skirt_path".
#
# The following table lists the configuration variables defined for each supported system, currently "cosma"
# and "obiwan".
#
#<TABLE>
#<TR><TD><B>Variable name</B></TD>  <TD><B>Description of value</B></TD></TR>
#<TR><TD>skirt_path</TD>            <TD>The absolute path to the SKIRT executable</TD></TR>
#<TR><TD>eagledata_path</TD>        <TD>A dictionary containing key-value pairs providing the absolute path to the
#                                       eagle data directory (value) for each relevant eagle simulation (key)</TD></TR>
#<TR><TD>catalogs_path</TD>         <TD>The absolute path to the directory containing the catalog files corresponding
#                                       to each of the the eagle snapshots used with this package</TD></TR>
#<TR><TD>database_path</TD>         <TD>The absolute path to the directory containing the SKIRT-runs database</TD></TR>
#<TR><TD>backup_path</TD>           <TD>The absolute path to the directory containing backups of the database</TD></TR>
#<TR><TD>templates_path</TD>        <TD>The absolute path to the directory containing ski file templates</TD></TR>
#<TR><TD>results_path</TD>          <TD>The absolute path to the directory containing the SKIRT results</TD></TR>
#<TR><TD>collections_path</TD>      <TD>The absolute path to the directory containing the collections with
#                                       statistics on sets of SKIRT simulation results</TD></TR>
#<TR><TD>plots_path</TD>            <TD>The absolute path to the directory containing the plots showing
#                                       statistics on sets of SKIRT simulation results</TD></TR>
#<TR><TD>default_eaglesim</TD>      <TD>The identifier of the eagle simulation currently in use (must be one of the
#                                       keys in the \em eagledata_path dictionary)</TD></TR>
#<TR><TD>default_redshift</TD>      <TD>The redshift of the snapshot currently in use</TD></TR>
#<TR><TD>queue</TD>                 <TD>The name of the queue to which jobs should be submitted, or None</TD></TR>
#<TR><TD>mpistyle</TD>              <TD>The style to invoke the mpirun command; one of 'generic' or 'lsf'</TD></TR>
#<TR><TD>nodes_per_job</TD>         <TD>The number of MPI computing nodes in each job (SKIRT-run);
#                                       specify 1 to disable MPI</TD></TR>
#<TR><TD>processes_per_node</TD>    <TD>The number of parallel MPI processes on each computing node;
#                                       specify 1 to disable MPI</TD></TR>
#<TR><TD>threads_per_process</TD>   <TD>The number of parallel threads in each MPI process;
#                                       specify zero to use the number of logical cores on the computing node</TD></TR>
#<TR><TD>maximum_hours</TD>         <TD>The maximum wall-time for each job (SKIRT-run) in hours;
#                                       this value is used (and enforced) only on batch queueing systems</TD></TR>
#</TABLE>
#
# The module also defines some variables identifying the current system, as listed in the following table.
#
#<TABLE>
#<TR><TD><B>Variable name</B></TD>  <TD><B>Description of value</B></TD></TR>
#<TR><TD>username</TD>              <TD>The login name of the user currently running the Python script</TD></TR>
#<TR><TD>hostname</TD>              <TD>The name of the host on which the Python script is running</TD></TR>
#</TABLE>
#
# Finally the module offers some utility functions for simple tasks such as obtaining a time stamp.
#

# -----------------------------------------------------------------

import datetime
import os
import os.path
import pwd
import socket

# -----------------------------------------------------------------

## This function returns a string representing the current time and date in the format "YYYY-MM-DD--hh-mm-ss".
# This format ensures proper collation (strings sort in date/time order), is easy to read for a human user,
# and can be used as part of a filename on any platform (since there are no nasty characters).
def timestamp():
    return datetime.datetime.now().strftime("%Y-%m-%d--%H-%M-%S")

## This function returns an absolute version of the specified path. The original path may be absolute, relative to
# a user's home folder, or relative to the current working directory.
def absolutepath(path):
    return os.path.realpath(os.path.expanduser(path))

# -----------------------------------------------------------------

# provide a dictionary of configuration dictionaries
configurations = {

    # the COSMA cluster in Durham
    'cosma': {
        'skirt_path': absolutepath("~/SKIRT/release/SKIRTmain/skirt"),
        'eagledata_path': { 'Ref100':  "/cosma5/data/Eagle/ScienceRuns/Planck1/L0100N1504/PE/REFERENCE/data/",
                           'Ref25':   "/cosma5/data/Eagle/ScienceRuns/Planck1/L0025N0752/PE/REFERENCE/data",
                           'Recal25': "/cosma5/data/Eagle/ScienceRuns/Planck1/L0025N0752/PE/RECALIBRATED/data" },
        'catalogs_path':  "/cosma5/data/Eagle/SkirtAnalysis/Catalogs",
        'database_path':  "/cosma5/data/Eagle/SkirtAnalysis/Database",
        'backup_path':    "/cosma5/data/Eagle/SkirtAnalysis/Backup",
        'templates_path': "/cosma5/data/Eagle/SkirtAnalysis/Templates",
        'results_path':   "/cosma5/data/Eagle/SkirtAnalysis/Results",
        'collections_path': "/cosma5/data/Eagle/SkirtAnalysis/Collections",
        'plots_path':     "/cosma5/data/Eagle/SkirtAnalysis/Plots",
        'default_eaglesim': 'Ref100',
        'default_redshift': 0,
        'queue': "cosma5",
        'mpistyle': 'lsf',
        'nodes_per_job': 8,
        'processes_per_node': 4,      # cosma5 nodes have 16 cores and 128GB of memory
        'threads_per_process': 4,
        'maximum_hours': 24,
    },

    # Peter's desktop at work
    'obiwan': {
        'skirt_path': absolutepath("~/SKIRT/release/SKIRTmain/skirt"),
        'eagledata_path': { 'Ref100': absolutepath("~/EAGLEsnapshots/L0100N1504REF"),
                           'Ref25': absolutepath("~/EAGLEsnapshots/L0025N0752REF"),
                           'Recal25': absolutepath("~/EAGLEsnapshots/L0025N0752RECAL") },
        'catalogs_path':  absolutepath("~/Dropbox/EAGLE/Catalogs"),
        'database_path':  absolutepath("~/Dropbox/EAGLE/Database"),
        'backup_path':    absolutepath("~/Dropbox/EAGLE/Backup"),
        'templates_path': absolutepath("~/Dropbox/EAGLE/Templates"),
        'results_path':   "/Volumes/EAGLEDATA/Results",
        'collections_path': absolutepath("~/Dropbox/EAGLE/Collections"),
        'plots_path':     absolutepath("~/Dropbox/EAGLE/Plots"),
        'default_eaglesim': 'Ref100',
        'default_redshift': 0,
        'queue': None,
        'mpistyle': 'generic',
        'nodes_per_job': 1,
        'processes_per_node': 2,
        'threads_per_process': 2,
        'maximum_hours': 0,
    }
}

# -----------------------------------------------------------------

# get the name of the user logged in on the terminal controlling this process; try various mechanisms
username = ""
if len(username)==0:
    try:
        username = os.environ['USER']
    except Exception:
        pass
if len(username)==0:
    try:
        username = os.getlogin()
    except Exception:
        pass
if len(username)==0:
    try:
        username = pwd.getpwuid(os.geteuid())[0]
    except Exception:
        pass

# get the name of the current host; try various mechanisms
hostname = ""
if len(hostname)==0:
    try:
        hostname = os.environ['LSB_QUEUE']      # for a queued job, use the queue name as host name
    except Exception:
        pass
if len(hostname)==0:
    try:
        hostname = os.environ['HOST']
    except Exception:
        pass
if len(hostname)==0:
    try:
        hostname = os.environ['HOSTNAME']
    except Exception:
        pass
if len(hostname)==0:
    try:
        if socket.gethostname().find('.')>=0:
            hostname=socket.gethostname()
        else:
            hostname=socket.gethostbyaddr(socket.gethostname())[0]
    except Exception:
        pass

# -----------------------------------------------------------------

# determine the configuration matching the current host
configuration = None
for configkey in configurations.keys():
    if configkey in hostname:
        configuration = configurations[configkey]
        break

if configuration == None:
    raise ValueError("Unknown host: " + hostname)

# -----------------------------------------------------------------

# copy all configurations to the module's namespace
globals().update(configurations)

# copy the settings for the current configuration to the module's namespace
globals().update(configuration)

# -----------------------------------------------------------------

# provide the table name prefix in the public EAGLE database corresponding to each of our shorthands
public_eagle_database_name = { 'Ref100':  "RefL0100N1504",
                               'Ref25':   "RefL0025N0752",
                               'Recal25': "RecalL0025N0752",
                             }

# get the user's account information for the public EAGLE database, if available
public_eagle_database_username = username
public_eagle_database_password = None
try:
    public_eagle_database_username, public_eagle_database_password = \
        open(os.path.join(database_path,"public_eagle_database_account_info.txt")).readline().split()
except Exception:
    pass

# -----------------------------------------------------------------
