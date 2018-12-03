
declare -A camera_from_site
camera_from_site+=( ["lsc"]="nres01" )
camera_from_site+=( ["elp"]="nres02" )
camera_from_site+=( ["cpt"]="nres03" )
camera_from_site+=( ["tlv"]="nres04" )

declare -A site_from_camera
site_from_camera+=( ["nres01"]="lsc" )
site_from_camera+=( ["nres02"]="elp" )
site_from_camera+=( ["nres03"]="cpt" )
site_from_camera+=( ["nres04"]="tlv" )
