  pro read_oci_geo_lut,lutfile,geo_lut

; Routine to read OCI geolocation LUT into a structure

;  Name		Type 	I/O	Description
;
;  lutfile	string	 I	Path/name of geolocation LUT file
;  geo_lut	structure O	Structure containing LUT parameters

; Open LUT file and read dimensions
  lutid = ncdf_open(lutfile)
  mdimid = ncdf_dimid(lutid,'number_of_HAM_sides')
  ncdf_diminq, lutid, mdimid, dimnm, msdim
  tdimid = ncdf_dimid(lutid,'number_of_tilts')
  ncdf_diminq, lutid, tdimid, dimnm, tdim
  bdimid = ncdf_dimid(lutid,'number_of_mce_boards')
  ncdf_diminq, lutid, bdimid, dimnm, bdim
  vdimid = ncdf_dimid(lutid,'vector_size')
  ncdf_diminq, lutid, vdimid, dimnm, vdim
  pdimid = ncdf_dimid(lutid,'planarity_coefficients')
  ncdf_diminq, lutid, pdimid, dimnm, pdim

; Define LUT structure
  geo_lut = { master_clock:0.d0, MCE_clock:0.d0, $
	sc_to_tilt:dblarr(vdim, vdim), $
	tilt_axis:dblarr(vdim), tilt_angles:dblarr(tdim), $
	tilt_home:0.d0, tilt_to_oci_mech:dblarr(vdim, vdim), $
	oci_mech_to_oci_opt:dblarr(vdim, vdim), $
	RTA_axis:dblarr(vdim), HAM_axis:dblarr(vdim), $
	HAM_AT_angles:dblarr(msdim), HAM_CT_angles:dblarr(msdim), $
	RTA_enc_scale:0.d0, HAM_enc_scale:0.d0, RTA_nadir:lonarr(bdim), $
	as_planarity:dblarr(pdim), at_planarity:dblarr(pdim) }

; Read parameters into structure
  tgid = ncdf_ncidinq(lutid,'time_params')
  ncdf_varget, tgid, 'master_clock', mac
  geo_lut.master_clock = mac
  ncdf_varget, tgid, 'MCE_clock', mcc
  geo_lut.MCE_clock = mcc

  cgid = ncdf_ncidinq(lutid,'coord_trans')
  ncdf_varget, cgid, 'sc_to_tilt', s2t
  geo_lut.sc_to_tilt(*,*) = s2t
  ncdf_varget, cgid, 'tilt_axis', tax 
  geo_lut.tilt_axis(*) = tax
  ncdf_varget, cgid, 'tilt_angles', tan
  geo_lut.tilt_angles(*) = tan
  ncdf_varget, cgid, 'tilt_home', thome
  geo_lut.tilt_home = thome
  ncdf_varget, cgid, 'tilt_to_oci_mech', t2o
  geo_lut.tilt_to_oci_mech(*,*) = t2o
  ncdf_varget, cgid, 'oci_mech_to_oci_opt', o2o
  geo_lut.oci_mech_to_oci_opt(*,*) = o2o

  rhgid = ncdf_ncidinq(lutid,'RTA_HAM_params')
  ncdf_varget, rhgid, 'RTA_axis', rax
  geo_lut.RTA_axis(*) = rax
  ncdf_varget, rhgid, 'HAM_axis', hax 
  geo_lut.HAM_axis(*) = hax
  ncdf_varget, rhgid, 'HAM_AT_angles', haa 
  geo_lut.HAM_AT_angles(*) = haa
  ncdf_varget, rhgid, 'HAM_CT_angles', hca 
  geo_lut.HAM_CT_angles(*) = hca
  ncdf_varget, rhgid, 'RTA_enc_scale', res
  geo_lut.RTA_enc_scale = res
  ncdf_varget, rhgid, 'HAM_enc_scale', hes
  geo_lut.HAM_enc_scale = hes
  ncdf_varget, rhgid, 'RTA_nadir', rtan
  geo_lut.RTA_nadir = rtan
  
  pgid = ncdf_ncidinq(lutid,'planarity')
  ncdf_varget, pgid, 'along_scan_planarity', asp
  geo_lut.as_planarity(*) = asp
  ncdf_varget, pgid, 'along_track_planarity', atp
  geo_lut.at_planarity(*) = atp


  ncdf_close,lutid

  return
  end

