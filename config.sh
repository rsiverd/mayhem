
arch_root="/archive/engineering"
save_root="/scratch/rsiverd/mayhem"
calib_conf="./calib_types.txt"
barrier_list="./barriers.txt"

## Camera storage locations:
#declare -A cam_storage
#cam_storage["nres"]="/archive/engineering/bpl/nres"
#cam_storage["nres01"]="/archive/engineering/lsc/nres01"

bias_stack_args="-r410 --itable"
dark_stack_args="-r410 --itable -c5"

#declare -A objects_keyword
#objects_keyword["thar01"]="thar&thar&none"
#objects_keyword["thar12"]="none&thar&thar"
#objects_keyword["tung01"]="tung&tung&none"
#objects_keyword["tung12"]="none&tung&tung"

##--------------------------------------------------------------------------##
## Calibration version requirements:
#min_biasvers="0.56"     # biases with version older than this are rebuilt
#min_darkvers="0.60"     #  darks with version older than this are rebuilt
#min_lampvers="0.35"     #  lamps with version older than this are rebuilt

## Image data version requirements for 'clean' images:
min_data_versions=( 0.58 0.63 0.42 )  # bias, dark, lamp

## Image code version requirements for stacked data:
min_code_versions=( 0.58 0.63 0.42 )  # bias, dark, lamp

######################################################################
# CHANGELOG (config.sh):
#---------------------------------------------------------------------
#
#  2018-08-03:
#     -- Increased min_data_versions and min_code_versions to force rebuild
#           of all old data now that nres-cdp-trim-oscan corrects gain.
#
#  2018-02-19:
#     -- After testing, removed commented-out code that had been moved.
#     -- Moved file/folder selection routines to new func/05_file_selection.sh.
#     -- Moved several date- and time-related functions to new source code
#           file func/03_time_and_date.sh. Functions migrated include:
#              date_is_listed(), dates_are_listed() to new
#
#  2018-01-05:
#     -- Introduced min_biasvers=0.55, min_darkvers=0.60, min_lampvers=0.35.
#     -- Moved barrier_check_pass() into separate aux/01_barriers.sh file.
#     -- update_output_header now also sets MVERSION for revision tracking.
#     -- Added this change log.
#

