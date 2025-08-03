declare -A all_options=(
  ['l2gen']="help= version= dump_options= dump_options_paramfile= dump_options_xmlfile= par= pversion= suite= eval= ifile= ilist= geofile= ofile= oformat= il2file= tgtfile= aerfile= metafile= l2prod= ctl_pt_incr= proc_ocean= proc_land= proc_sst= proc_cloud= atmocor= mode= aer_opt= aer_wave_short= aer_wave_long= aer_wave_base= aer_swir_short= aer_swir_long= aer_rrs_short= aer_rrs_long= aermodmin= aermodmax= aermodrat= aer_angstrom= aer_iter_max= mumm_alpha= mumm_gamma= mumm_epsilon= absaer_opt= glint_opt= oxaband_opt= cirrus_opt= filter_opt= filter_file= brdf_opt= fqfile= parfile= gas_opt= atrem_opt= atrem_full= atrem_geom= atrem_model= atrem_splitpaths= iop_opt= cphyt_opt= seawater_opt= polfile= pol_opt= band_shift_opt= giop_aph_opt= giop_aph_file= giop_aph_s= giop_adg_opt= giop_adg_file= giop_adg_s= giop_bbp_opt= giop_bbp_file= giop_bbp_s= giop_acdom_opt= giop_acdom_file= giop_anap_opt= giop_anap_file= giop_bbph_opt= giop_bbph_file= giop_bbnap_opt= giop_bbnap_file= giop_rrs_opt= giop_rrs_diff= giop_grd= giop_wave= giop_maxiter= giop_fit_opt= gsm_opt= gsm_fit= gsm_adg_s= gsm_bbp_s= gsm_aphw= gsm_aphs= qaa_adg_s= qaa_wave= chloc2_wave= chloc2_coef= chloc3_wave= chloc3_coef= chloc4_wave= chloc4_coef= kd2_wave= kd2_coef= flh_offset= vcnnfile= picfile= owtfile= owtchlerrfile= aermodfile= uncertaintyfile= aermodels= met1= met2= met3= ozone1= ozone2= ozone3= rad1= rad2= rad3= anc_profile1= anc_profile2= anc_profile3= anc_aerosol1= anc_aerosol2= anc_aerosol3= sfc_albedo= anc_cor_file= pixel_anc_file= land= water= demfile= dem_auxfile= mldfile= icefile= ice_threshold= sstcoeffile= dsdicoeffile= sstssesfile= sst4coeffile= sst4ssesfile= sst3coeffile= sst3ssesfile= sstfile= sstreftype= sstrefdif= viirsnv7= viirsnosisaf= no2file= alphafile= tauafile= flaguse= xcalbox= xcalboxcenter= xcalpervalid= xcalsubsmpl= chlthreshold= aotthreshold= coccolith= cirrus_thresh= taua= absaer= rhoamin= epsmin= epsmax= tauamax= nLwmin= wsmax= windspeed= windangle= pressure= ozone= relhumid= watervapor= vcal_opt= vcal_chl= vcal_solz= vcal_nLw= vcal_Lw= vcal_depth= vcal_min_nbin= vcal_min_nscene= owmcfile= north= south= east= west= xbox= ybox= subsamp= prodxmlfile= breflectfile= bpar_validate_opt= bpar_elev_opt= bpar_elev_value= deflate= raman_opt= gmpfile= water_spectra_file= shallow_water_depth= avw_coef= cloud_hgt_file= doi= wavelength_3d= mbac_wave= calfile= rad_opt= viirscalparfile= geom_per_band= xcalfile= xcal_opt= xcal_wave= resolution= newavhrrcal= ch22detcor= ch23detcor= sl_pixl= sl_frac= outband_opt= maskland= maskbath= maskcloud= maskglint= masksunzen= masksatzen= maskhilt= maskstlight= sunzen= satzen= hipol= glint_thresh= cloud_thresh= cloud_wave= cloud_eps= cloud_mask_file= offset= gain= spixl= epixl= dpixl= sline= eline= dline="
  ['l2bin']="help= version= dump_options= dump_options_paramfile= dump_options_xmlfile= ifile= par= ofile= suite= eval= resolution= resolve= deflate= verbose= latnorth= latsouth= loneast= lonwest= fileuse= delta_crossing_time= qual_prod= sday= eday= night= flaguse= l3bprod= area_weighting= output_wavelengths= oprodname= composite_scheme= composite_prod= pversion= prodtype= qual_max= rowgroup= minobs="
  ['l3mapgen']="help= version= dump_options= dump_options_paramfile= dump_options_xmlfile= par= suite= ifile= ofile= oformat= ofile_product_tag= ofile2= oformat2= deflate= product= wavelength_3d= resolution= width= projection= write_projtext= central_meridian= lat_ts= lat_0= lat_1= lat_2= azimuth= utm_zone= north= south= east= west= trimNSEW= interp= apply_pal= palfile= use_transparency= datamin= datamax= scale_type= quiet= pversion= use_quality= quality_product= use_rgb= product_rgb= fudge= threshold= num_cache= mask_land= rgb_land= land= full_latlon= doi=" 
  ['l1bgen_oci']="help= version= dump_options= dump_options_paramfile= dump_options_xmlfile= par= ifile= ofile= cal_lut= geo_lut= doi= pversion= demfile= radiance= disable_geo= ephfile="
  ['focs']="process include="
  ['l1cgen']="help= version= dump_options= dump_options_paramfile= dump_options_xmlfile= par= ifile= l1c_grid= l1c_anc= ofile= outlist= verbose= pversion= doi= l2prod= mode= south= north= west= east= selgran= ix_l1cprod= selday= selmon= selyear= grid_resolution= sensor= gransize= grantype= bintype= start_time= end_time= history= start_timeflag= swath_num= projection= sort_method= demcloud_flag= cloud_height= terrain_correct= cloud_correct= cloud_type= overlap_vflag= overlap_pflag= overlap_bflag= unc_meth= unc_thres_v= unc_thres_p= unc_thres_b= demfile="
  ['l3bindump']="help= version= dump_options= dump_options_paramfile= dump_options_xmlfile= par= ifile= ofile= oformat= l3bprod= bin_number= north= south= east= west= lat= lon= radius= verbose="
  ["l3mapmerge"]="ifile ofile --version --product --doi --pversion"
  ['l3bin']="help= version= verbose= dump_options= dump_options_paramfile= dump_options_xmlfile= par= pversion= ifile= ofile= oformat= merged= latnorth= latsouth= loneast= lonwest= sday= eday= deflate= orbit1= orbit2= median= noext= unit_wgt= composite_scheme= composite_prod= reduce_fac= resolve= prod= doi="
  ['l2extract']="help= version= dump_options= dump_options_paramfile= dump_options_xmlfile= par= ifile= ofile= product= spix= epix= sline= eline= suite= verbose= wavelist="
  ['get_product_info']="help= l= r= csv= version= dump_options= dump_options_paramfile= dump_options_xmlfile= sensor="
  ['vcalmerge']="help= version= dump_options= dump_options_paramfile= dump_options_xmlfile= par= pversion= ifile= ifile1= ifile2= ofile= ofile1= ofile2= oformat= spixl= epixl= sline= eline= flaguse= deflate="
  ['l2binmatch']="help= version= dump_options= dump_options_paramfile= dump_options_xmlfile= par= ifile= ofile= l2prod= flaguse= oformat= deflate= tgtfile= spixl= epixl= dpixl= sline= eline= dline= subsamp= band_shift_opt= vcal_depth= vcal_min_nbin= vcal_min_nscene= demfile="
  ['ancgen']="help= version= dump_options= dump_options_paramfile= dump_options_xmlfile= par= targetfile= ancfile= merraprofile= merrametfile= merraaerfile= geoscffile= cldmask1= cldmask2= cldmask3= cldprod1= cldprod2= cldprod3= albedofile= chlfile= ch4file= co2file= n2ofile= gebcofile="
  ['l3binmerge']="help= version= verbose= dump_options= dump_options_paramfile= dump_options_xmlfile= par= pversion= ifile= ofile= latnorth= latsouth= loneast= lonwest= noext= union_bins= prod="
  ['l1det2det']="help= version= dump_options= dump_options_paramfile= dump_options_xmlfile= par= pversion= suite= ifile= ofile= oformat= geofile= l2prod= ybox= chlthreshold= aotthreshold= cloud_thresh= flaguse= glint_thresh= spixl= epixl= sline= eline="
  ['l1bgeneric']="version= dump_options= dump_options_paramfile= dump_options_xmlfile= par= pversion= suite= ifile= ofile= oformat= fqfile= parfile= georegion_file= calfile= xcalfile= sl_pixl= sl_frac= gain= offset= spixl= epixl= sline= eline="
  ['l1agen_oci']="help= version= dump_options= dump_options_paramfile= dump_options_xmlfile= par= ifile= ofile= granule_len= start_time= maxgap= hktlist= swir_loff_set= outlist= nametag= doi= pversion= noSPW="
  ['val_extract']="ofile= dump= dump_type= ifile= sline= spixl= boxsize= boxsizekm= slat= slon= elat= elon= global_att= variable_att= l2qc= valid_ranges= ignore_flags= skip_stats= count_flags="
)

declare -A python_scripts=(
    ["modis_L1B"]="-h --version -p --parfile -o --okm -k --hkm -q --qkm -c --obc -l --lutver -d --lutdir -x --del-okm -y --del-hkm -z --del-qkm --keep-obc -v --verbose --log"
    ["find_matchup"]="-h --sat --data_type --slat --elat --slon --elon --stime --max_time_diff --etime --seabass_file --get_data --verbose"
    ["l1info"]="-d -i -n -s -v -o"
    ["l3mapmerge"]="ifile ofile --version --product --doi --pversion"
    ["nccmp"]="-A --Attribute -a --all -e --extent -b --verbose -d --data -C --maxdiff -D --debug -f --force -F --fortran -g --global -G --globalex -h --history -H --help --usage -m --metadata -M --missing -N --nans-are-equal -n --notolerance -p --precision -s --report-identical-files -t --tolerance -T --Tolerance -v --variable -V --version -x --exclude -w --warn"
)

# Define valid values for specific options
declare -A option_values=(
  # Existing options
  ['area_weighting']="0 1 2 3"
  ['suite']="SST NSST IOP OC KD RRS RRS_SB CLD CLDMASK PAR"
  ['resolve']="H Q HQ HH 1 2 4 9 18 36 1D HD QD"
  ['resolution']="0.5km 250m 100m 50m 1.1km 2.3km 4.6km 9.2km 18.5km 36km 1 degree 0.5 degree 0.25 degree"
  ['projection']="smi platecarree mollweide lambert albersconic mercator transmerc utm obliquemerc ease2 ortho stere conus alaska gibs raw"
  ['eval']="0 1 2 16 32 64 128 256 1024 2048 4096 8192 32768"
  ['mode']="0 1 2 3"
  ['aer_opt']="-99 0 -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -17 -18 -19"  # Simplified, assuming >0 is free-form model number
  ['glint_opt']="0 1 2"
  ['absaer_opt']="0 1 2"
  ['oxaband_opt']="0 1 2"
  ['proc_ocean']="0 1 2"  # Updated to include 2 from help
  ['proc_land']="0 1"
  ['proc_sst']="0 1"
  ['cirrus_opt']="0 1"
  ['filter_opt']="0 1"

  # New l2gen options from help menu
  ['help']="0 1"
  ['verbose']="0 1"
  ['dump_options']="0 1"
  ['oformat']="netCDF4 hdf4"
  ['proc_cloud']="0 1"
  ['atmocor']="0 1"
  ['brdf_opt']="0 1 3 7 15 19"
  ['gas_opt']="0 1 2 4 8 16 32 64 128 256"
  ['atrem_opt']="0 1 2 4 8 16 32 64"
  ['atrem_full']="0 1"
  ['atrem_geom']="0 1"
  ['atrem_model']="0 1 2 3 4 5 6"
  ['atrem_splitpaths']="0 1"
  ['iop_opt']="0 1 2 3 4 5 6 7"
  ['cphyt_opt']="1 2"
  ['seawater_opt']="0 1"
  ['pol_opt']="-1 0 1 2 3 4"
  ['band_shift_opt']="0 1"
  ['giop_aph_opt']="0 2 3"
  ['giop_adg_opt']="0 1 2 3"
  ['giop_bbp_opt']="0 1 2 3 5 6 7 8 9 10"
  ['giop_acdom_opt']="0 1"
  ['giop_anap_opt']="0 1"
  ['giop_bbph_opt']="0 1"
  ['giop_bbnap_opt']="0 1"
  ['giop_rrs_opt']="0 1"
  ['giop_fit_opt']="0 1 3 4"
  ['gsm_opt']="0 1"
  ['gsm_fit']="0 1"
  ['sstreftype']="0 1 2 3 4 5 6 7 8"
  ['viirsnv7']="-1 1"
  ['viirsnosisaf']="0 1"
  ['bpar_validate_opt']="0 1"
  ['bpar_elev_opt']="0 1"
  ['deflate']="0 1 2 3 4 5 6 7 8 9"  # Assuming standard deflate levels 0-9
  ['raman_opt']="0 1 2 3"
  ['rad_opt']="0 1"
  ['geom_per_band']="0 1"
  ['resolution']="-1 250 500 1000"  # l2gen-specific resolution values
  ['newavhrrcal']="0 1"
  ['outband_opt']="0 2 99"
  ['maskland']="0 1"
  ['maskbath']="0 1"
  ['maskcloud']="0 1"
  ['maskglint']="0 1"
  ['masksunzen']="0 1"
  ['masksatzen']="0 1"
  ['maskhilt']="0 1"
  ['maskstlight']="0 1"

    # New l3mapgen options
  ['oformat']="netcdf4 hdf4 png ppm tiff"
  ['oformat2']="netcdf4 hdf4 png ppm tiff"
  ['deflate']="0 1 2 3 4 5 6 7 8 9"  # Standard deflate levels
  ['product']="chlor_a Kd_490:avg chlor_a:stdev chlor_a:nobs chlor_a:nscenes chlor_a:obs_time chlor_a:bin_num"  # Example products with modifiers
  ['interp']="nearest bin area"
  ['apply_pal']="0 1"
  ['use_transparency']="0 1"
  ['scale_type']="linear log arctan"
  ['quiet']="0 1"
  ['use_quality']="0 1"
  ['use_rgb']="0 1"
  ['write_projtext']="0 1"
  ['trimNSEW']="0 1"
  ['mask_land']="0 1"
  ['full_latlon']="0 1"
  # New l1bgen_oci options
  ['radiance']="0 1"
  ['disable_geo']="0 1"
  # New l3bin options
  ['median']="0 1"  # Assuming boolean-like behavior (0=off, 1=on)
  ['noext']="0 1"
  ['oformat']="netcdf4 hdf4 hdf5"
  ['unit_wgt']="0 1"  # Assuming boolean-like behavior (0=off, 1=on)
  ['composite_scheme']="min max"
  ['reduce_fac']="1 2 4 8 16 32 64 128 256"  # Power of 2 values
# New l1cgen options
  ['mode']="5 7 8"  # Specific to l1cgen
  ['ix_l1cprod']="0 1 2"  # 0: pc, 1: vsf, 2: dpr
  ['sensor']="SPEXONE OCI HARP2 MISR"
  ['grantype']="0 1"  # 0: granules, 1: swath
  ['bintype']="0 1"  # 0: discrete, 1: area-weighting
  ['start_timeflag']="0 1"  # 0: time zero, 1: starting coverage
  ['swath_num']="1 2"  # 1: ascending, 2: descending
  ['projection']="0 1"  # 0: SOCEA, 1: SOCEA-2 (overrides l3mapgen values for l1cgen)
  ['sort_method']="0 1"  # 0: Orbital-vectors, 1: SADDLEBACK SEARCH
  ['demcloud_flag']="0 1"  # 0: geoid/L1C height, 1: orthometric/dem height
  ['terrain_correct']="0 1"  # 0: off, 1: on
  ['cloud_correct']="0 1 2"  # 0: no correction, 1: constant CTH from ANC, 2: constant CTH from L1C grid
  ['cloud_type']="0 1"  # 0: water, 1: ice
  ['overlap_vflag']="0 1"  # 0: off, 1: on
  ['overlap_pflag']="0 1"  # 0: off, 1: on
  ['overlap_bflag']="0 1"  # 0: off, 1: on
  ['unc_meth']="0 1"  # 0: error propagation, 1: Monte Carlo
  ['band_shift_opt']="0 1"  # 0: linear interpolation, 1: bio-optical bandshift
  ['noext']="0 1"
  ['union_bins']="0 1"
  ['noSPW']="0 1"
  ['start_time']="YYYYmmddTHHMMSS YYYY-mm-ddTHH:MM:SS"
  ['swir_loff_set']="0,0,0,0,0,0,0,0,0"  # Example format
)
# program specific options
# modes
declare -A modes_options=(
    ["l2gen"]="0 1 2 3"
    ["l1cgen"]="5 7 8" 
)

# formats
declare -A oformat_options=(
    ["l2gen"]="netCDF4 hdf4"
    ["l3mapgen"]="netcdf4 hdf4 png ppm tiff"
    ["l3bin"]="netcdf4 hdf4 hdf5"
    ["vcalmerge"]="netcdf4 hdf4"
    ["l2binmatch"]="netcdf4 hdf4"
)

# declare
declare -A resolution_options=(
    ["l2gen"]="-1 250 500 1000"
    ["l3mapgen"]="90km 36km 18km 9km 4km 2km 1km hkm qkm smi smi4 land"  
    ["l2bin"]="0.5km 250m 100m 50m 1.1km 2.3km 4.6km 9.2km 18.5km 36km 1 degree 0.5 degree 0.25 degree"
)


# declare
declare -A projection_options=(
     ["l3mapgen"]="smi platecarree mollweide lambert albersconic mercator transmerc utm obliquemerc ease2 ortho stere conus alaska gibs raw" #"0.5km 250m 100m 50m 1.1km 2.3km 4.6km 9.2km 18.5km 36km 1 degree 0.5 degree 0.25 degree"
     ["l1cgen"]="0 1"
)
_default_completion() {
    local cur prev opts

    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    opts=${all_options[$1]}

    # Explicitly handle area_weighting= to suggest 0-3
    if [[ -n "${option_values[$prev]+x}"  && "$cur" == "=" ]]; then
        if [[  "$prev" == "mode" ]]; then
            COMPREPLY=( $(compgen -W "${modes_options[$1]}") )
            return 0
        fi
        if [[  "$prev" == "oformat" ]]; then
            COMPREPLY=( $(compgen -W "${oformat_options[$1]}") )
            return 0
        fi
        if [[  "$prev" == "resolution" ]]; then
            COMPREPLY=( $(compgen -W "${resolution_options[$1]}") )
            return 0
        fi
        if [[  "$prev" == "projection" ]]; then
            COMPREPLY=( $(compgen -W "${projection_options[$1]}") )
            return 0
        fi          
        COMPREPLY=( $(compgen -W "${option_values[$prev]}") )
        return 0
    fi

    # General options completion
    COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )

    # Filename completion after '='
    if [[ "${prev}" == '=' ]]; then
        cur="${cur//\\ / }"
        [[ "${cur}" == "~/"* ]] && cur="${cur/\~/$HOME}"
        compopt -o filenames
        compopt +o nospace
        local files=("${cur}"*)
        if [[ -e "${files[0]}" ]]; then
            COMPREPLY=( "${files[@]// /\ }" )
        else
            COMPREPLY=()
        fi
        return 0
    fi

    # Trigger filename completion when cursor is on '='
    if [[ "${cur}" == '=' ]]; then
        compopt -o default
        compopt +o nospace
        local files=(*)
        if [[ -e "${files[0]}" ]]; then
            COMPREPLY=( "${files[@]// /\ }" )
        else
            COMPREPLY=()
        fi
        return 0
    fi
}

_program_autocomplete() {
    local cur prev opts
    COMPREPLY=()  # Array for completion results
    cur="${COMP_WORDS[COMP_CWORD]}"  # Current word being completed
    prev="${COMP_WORDS[COMP_CWORD-1]}"  # Previous word
    
    # Define the available options
    opts=${python_scripts[$1]}
    
    # Generate completion matches
    COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
    
    return 0
}

for i in "${!all_options[@]}"
do
  eval "_completion_${i}() { _default_completion ${i}; }"
  complete -o nospace -F "_completion_${i}" "${i}"
done

for i in "${!python_scripts[@]}"
do
  eval "_completion_${i}() { _program_autocomplete ${i}; }"
  complete -F "_completion_${i}" "${i}"
done

_install_ocssw_autocomplete() {
    local cur prev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # All available options
    opts="-h --help --version --list_tags --installed_tag --status --update --diff --difftool \
          --create_change --deliver_change --deliver_email_from -t --tag -i --install_dir \
          -b --base_url -l --local_dir -s --save_dir -c --clean --wget -a --arch -v --verbose \
          --bin --opt --src --luts --viirs_l1_bin --root --python --opt_src --afrt --aquaverse \
          --avhrr --aviris --common --czcs --eval --goci --harp2 --hawkeye --hico --l5tm \
          --l7etmp --meris --misr --modisa --modist --msis2a --msis2b --oci --ocm1 --ocm2 \
          --ocrvc --octs --olcis3a --olcis3b --olil8 --olil9 --prism --sabiamar --seawifs \
          --sgli --spexone --viirsdem --viirsj1 --viirsj2 --viirsn --wv3 --aerosol --cloud \
          --telemetery --benchmark --viirs_l1_benchmark --direct_broadcast --seadas --odps \
          --viirs_l1 --all"

    # Handle options that take arguments
    case "$prev" in
        --difftool|-i|--install_dir|-b|--base_url|-l|--local_dir|-s|--save_dir|--deliver_email_from)
            # For directory options, suggest files/directories
            if [[ "$prev" =~ "(--install_dir|-i|--local_dir|-l|--save_dir|-s|--difftool)" ]]; then
                COMPREPLY=( $(compgen -f -- "$cur") )
            fi
            return 0
            ;;
        -a|--arch)
            # Specific architecture options
            COMPREPLY=( $(compgen -W "linux_64 linux_hpc macosx_intel macosx_arm64 odps" -- "$cur") )
            return 0
            ;;
        -t|--tag)
            # look up the avaliable tags
            COMPREPLY=( $(compgen -W "`install_ocssw --list_tags`" -- "$cur") )
            return 0
            ;;
    esac

    # Default case: suggest all options
    COMPREPLY=( $(compgen -W "$opts" -- "$cur") )
    return 0
}

# Register the completion
complete -F _install_ocssw_autocomplete install_ocssw


_obdaac_download_autocomplete() {
    local cur prev opts
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # Debug output to see what's happening
    # echo "DEBUG: cur='$cur' prev='$prev' COMP_CWORD='$COMP_CWORD'" >&2

    # Minimal set of options
    opts="-h --help -v --verbose --filelist --odir"

    # Default: suggest all options
    COMPREPLY=( $(compgen -W "$opts" -- "$cur") )
    return 0
}

# Register completion
complete -F _obdaac_download_autocomplete obdaac_download
