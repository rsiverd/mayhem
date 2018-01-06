## NOTE: this command-line parsing code is SOURCED by the caller because it
## needs to run in the script's global namespace.

##--------------------------------------------------------------------------##
## Parse command line with getopt (reorders and stores CL args):
s_opts="A:B:cko:rhqtv"
l_opts="after:,before:,clobber,keep_clean,output:,random"
l_opts+=",thar01,thar12,tung01,tung12"
l_opts+=",debug,deps,help,quiet,timer,verbose"
args=`getopt -n $this_prog -o $s_opts -l $l_opts -- "$@"` ; failed=$?

## Check for parsing errors:
if [ $failed -ne 0 ]; then
   echo "Try \`$this_prog --help' for more information." >&2
   exit 1
fi

## Change the arguments in the environment (e.g., $1, $2, ... $N):
eval set -- "$args"

## Loop through arguments (shift removes them from environment):
while true ; do
   case $1 in
      #-------------------------------------------------------------------
     #--START)
     #   confirmed=1
     #   shift
     #   ;;
      #-------------------------------------------------------------------
      -A|--after)
         case $2 in
            -*)
               msg="Option -A|--after needs a non-negative integer argument!"
               Recho "\n${msg}\n" >&2
               usage >&2
               exit 1
               ;;
            *)
               if !( int_check_pass $2 ) || [ $2 -lt 0 ]; then
                  Recho "Invalid days_after: " >&2 ; Yecho "$2 \n\n" >&2
                  exit 1
               fi
               days_after=$2
               ;;
         esac
         [ $vlevel -ge 0 ] && yecho "Stack days after: ${days_after}\n" >&2
         shift 2
         ;;
      #-------------------------------------------------------------------
      -B|--before)
         case $2 in
            -*)
               msg="Option -B|--before needs a non-negative integer argument!"
               Recho "\n${msg}\n" >&2
               usage >&2
               exit 1
               ;;
            *)
               if !( int_check_pass $2 ) || [ $2 -lt 0 ]; then
                  Recho "Invalid days_prior: " >&2 ; Yecho "$2 \n\n" >&2
                  exit 1
               fi
               days_prior=$2
               ;;
         esac
         [ $vlevel -ge 0 ] && yecho "Stack days before: ${days_prior}\n" >&2
         shift 2
         ;;
      #-------------------------------------------------------------------
      -c|--clobber)
         [ $vlevel -ge 0 ] && yecho "Enabled output clobber!\n" >&2
         clobber=1
         shift
         ;;
      #-------------------------------------------------------------------
      -k|--keep_clean)
         [ $vlevel -ge 0 ] && yecho "Preserving 'clean' files!\n" >&2
         keep_clean=1
         shift
         ;;
      #-------------------------------------------------------------------
     #-n|--number)
     #   case $2 in
     #      -*)
     #         msg="Option -n|--number requires an argument!"
     #         #msg="Option -n|--number needs a positive integer argument!"
     #         #msg="Option -n|--number needs a positive numerical argument!"
     #         Recho "\n${msg}\n" >&2
     #         usage >&2
     #         exit 1
     #         ;;
     #      *)
     ###       if !( num_check_pass $2 ); then
     ###       if !( num_check_pass $2 ) || (is_negative $2); then
     ###       if !( int_check_pass $2 ) || [ $2 -lt 0 ]; then
     #            Recho "Invalid value: " >&2 ; Yecho "$2 \n\n" >&2
     #            exit 1
     #         fi
     #         num_val=$2
     #         ;;
     #   esac
     #   [ $vlevel -ge 0 ] && yecho "Using value: ${num_val}\n" >&2
     #   shift 2
     #   ;;
      #-------------------------------------------------------------------
      -o|--output)
         case $2 in
            -*)
               Recho "\nOption -o|--output requires an argument!\n" >&2
               usage >&2
               exit 1
               ;;
            *)
               save_file=$2
               # check value here ...
               if [ $clobber -eq 0 ] && [ -f $save_file ]; then
                  Recho "\nFile already exists: " >&2
                  Yecho "$save_file \n\n" >&2
                  exit 1
               fi
               ;;
         esac
         [ $vlevel -ge 0 ] && yecho "Output to: $save_file \n" >&2
         shift 2
         ;;
      -r|--random)
         [ $vlevel -ge 0 ] && yecho "Randomizing order!\n" >&2
         shuffle=1
         shift
         ;;
      #-------------------------------------------------------------------
      # Lamp type and channel selection:
      --thar01|--thar12|--tung01|--tung12)
         calib_type="${1##--}"
         [ $vlevel -ge 0 ] && yecho "Calib type: ${calib_type}\n" >&2
         shift
         ;;
      #-------------------------------------------------------------------
      # Additional options (output control etc.):
      --debug)
         yecho "Debugging mode enabled!\n"
         debug=1
         shift
         ;;
      --deps)
         echo ${need_exec[*]}
         exit 0
         ;;
     #-f|--force)
     #   [ $vlevel -ge 0 ] && yecho "Output clobbering enabled!\n" >&2
     #   [ $vlevel -ge 0 ] && yecho "Forcing <WHAT THIS SCRIPT DOES>!\n" >&2
     #   clobber=1
     #   force=1
     #   shift
     #   ;;
      -h|--help)
         usage
         exit 0
         ;;
      -q|--quiet)
         (( vlevel-- ))
         shift
         ;;
      -t|--timer)
         [ $vlevel -ge 1 ] && yecho "Timing script execution!\n" >&2
         timer=1
         shift
         ;;
      -v|--verbose)
         (( vlevel++ ))
         shift
         ;;
      #-------------------------------------------------------------------
      --)
         shift
         break
         ;;
      *)
         echo -e "\n$this_prog error: unhandled option: $1 \n" >&2
         exit 1
         ;;
      #-------------------------------------------------------------------
   esac
done

