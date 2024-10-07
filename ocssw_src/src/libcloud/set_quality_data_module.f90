 module set_quality_data_module
 
 implicit none
 
 contains
	 
 subroutine set_quality_data(xsize, ysize)

   use core_arrays
   use libraryarrays, only: water_radii, ice_radii, number_waterradii, number_iceradii
   use ch_xfr, only: c2_sensor_id, OCI_ID, OCIS_ID

   integer , intent(in) :: xsize, ysize

   integer :: i, j

   integer, parameter :: marginal = 1, good = 2, very_good = 3, no_confidence = 0
   
   integer*1, parameter :: no_cloud_mask = 0, no_cloud = 1, liquid = 2, &
   							ice = 3, unknown = 4, good_retrieval = 8

   do i = 1, xsize
     do j = 1, ysize

       if (cloudsummary(i,j)%cloudmask_determined) then

            processing_information(i,j)%path_and_outcome = no_cloud
            processing_information(i,j)%path_and_outcome_PCL = no_cloud
            processing_information(i,j)%path_and_outcome_16 = no_cloud
            processing_information(i,j)%path_and_outcome_16_PCL = no_cloud
            processing_information(i,j)%path_and_outcome_37 = no_cloud
            processing_information(i,j)%path_and_outcome_37_PCL = no_cloud
            processing_information(i,j)%path_and_outcome_22 = no_cloud
            processing_information(i,j)%path_and_outcome_22_PCL = no_cloud
       
          if(cloudsummary(i,j)%watercloud) then
            processing_information(i,j)%path_and_outcome = liquid
            processing_information(i,j)%path_and_outcome_PCL = liquid
            processing_information(i,j)%path_and_outcome_16 = liquid
            processing_information(i,j)%path_and_outcome_16_PCL = liquid
            processing_information(i,j)%path_and_outcome_37 = liquid
            processing_information(i,j)%path_and_outcome_37_PCL = liquid
            processing_information(i,j)%path_and_outcome_22 = liquid
            processing_information(i,j)%path_and_outcome_22_PCL = liquid
          elseif(cloudsummary(i,j)%icecloud) then
            processing_information(i,j)%path_and_outcome = ice
            processing_information(i,j)%path_and_outcome_PCL = ice
            processing_information(i,j)%path_and_outcome_16 = ice
            processing_information(i,j)%path_and_outcome_16_PCL = ice
            processing_information(i,j)%path_and_outcome_37 = ice
            processing_information(i,j)%path_and_outcome_37_PCL = ice
            processing_information(i,j)%path_and_outcome_22 = ice
            processing_information(i,j)%path_and_outcome_22_PCL = ice
          elseif(cloudsummary(i,j)%unknowncloud) then
            processing_information(i,j)%path_and_outcome = unknown
            processing_information(i,j)%path_and_outcome_PCL = unknown
            processing_information(i,j)%path_and_outcome_16 = unknown
            processing_information(i,j)%path_and_outcome_16_PCL = unknown
            processing_information(i,j)%path_and_outcome_37 = unknown
            processing_information(i,j)%path_and_outcome_37_PCL = unknown
            processing_information(i,j)%path_and_outcome_22 = unknown
            processing_information(i,j)%path_and_outcome_22_PCL = unknown
          endif
       else 
          processing_information(i,j)%path_and_outcome = no_cloud_mask
          processing_information(i,j)%path_and_outcome_PCL = no_cloud_mask
          processing_information(i,j)%path_and_outcome_16 = no_cloud_mask
          processing_information(i,j)%path_and_outcome_16_PCL = no_cloud_mask
          processing_information(i,j)%path_and_outcome_37 = no_cloud_mask
          processing_information(i,j)%path_and_outcome_37_PCL = no_cloud_mask
          processing_information(i,j)%path_and_outcome_22 = no_cloud_mask
          processing_information(i,j)%path_and_outcome_22_PCL = no_cloud_mask

       endif 

       if (cloudsummary(i,j)%cloudmask_determined .and. &
            ( cloudsummary(i,j)%ocean_surface .or.      &
              cloudsummary(i,j)%snowice_surface)      )  then

            processing_information(i,j)%path_and_outcome_1621 = no_cloud
            processing_information(i,j)%path_and_outcome_1621_PCL = no_cloud
			
          if(cloudsummary(i,j)%watercloud) then
            processing_information(i,j)%path_and_outcome_1621 = liquid
            processing_information(i,j)%path_and_outcome_1621_PCL = liquid
          elseif(cloudsummary(i,j)%icecloud) then
            processing_information(i,j)%path_and_outcome_1621 = ice
            processing_information(i,j)%path_and_outcome_1621_PCL = ice
          elseif(cloudsummary(i,j)%unknowncloud) then
            processing_information(i,j)%path_and_outcome_1621 = unknown
            processing_information(i,j)%path_and_outcome_1621_PCL = unknown
          endif 
		  
       else
          processing_information(i,j)%path_and_outcome_1621 = no_cloud_mask
          processing_information(i,j)%path_and_outcome_1621_PCL = no_cloud_mask
       endif

#ifdef VIIRS_INST 
      if (processing_information(i,j)%spectral_VNSWIR_21 == 0) &
      	processing_information(i,j)%path_and_outcome = no_cloud_mask
#endif      
      
		!  confidence = 11, usefulness = 1
       	if (optical_thickness_final(i,j) > 0. .and. effective_radius_21_final(i,j) > 0.) then 
       		processing_information(i,j)%optical_thickness_GC = 7
       		processing_information(i,j)%effective_radius_GC = 7
       		processing_information(i,j)%water_path_GC = 7
       	
       	    processing_information(i,j)%path_and_outcome = &
         				processing_information(i,j)%path_and_outcome + good_retrieval  
	
       	else 
       		processing_information(i,j)%optical_thickness_GC = 0
       		processing_information(i,j)%effective_radius_GC = 0
       		processing_information(i,j)%water_path_GC = 0
		endif       		

		!  confidence = 11, usefulness = 1
       	if (optical_thickness_final_PCL(i,j) > 0. .and. effective_radius_21_final_PCL(i,j) > 0.) then 
       		processing_information(i,j)%optical_thickness_GC = 7
       		processing_information(i,j)%effective_radius_GC = 7
       		processing_information(i,j)%water_path_GC = 7

         	processing_information(i,j)%path_and_outcome_PCL = &
         				processing_information(i,j)%path_and_outcome_PCL + good_retrieval  
       	endif 
     		

		!  confidence = 11, usefulness = 1
       	if (optical_thickness_1621_final(i,j) > 0. .and. effective_radius_1621_final(i,j) > 0.)  then
      		processing_information(i,j)%optical_thickness_1621_GC = 7
       		processing_information(i,j)%effective_radius_1621_GC = 7
       		processing_information(i,j)%water_path_1621_GC = 7

	        processing_information(i,j)%path_and_outcome_1621 = &
     					processing_information(i,j)%path_and_outcome_1621 + good_retrieval  
		else 
      		processing_information(i,j)%optical_thickness_1621_GC = 0
       		processing_information(i,j)%effective_radius_1621_GC = 0
       		processing_information(i,j)%water_path_1621_GC = 0
		endif

		!  confidence = 11, usefulness = 1
       	if (optical_thickness_1621_final_PCL(i,j) > 0. .and. effective_radius_1621_final_PCL(i,j) > 0.)  then
      		processing_information(i,j)%optical_thickness_1621_GC = 7
       		processing_information(i,j)%effective_radius_1621_GC = 7
       		processing_information(i,j)%water_path_1621_GC = 7
         	
         	processing_information(i,j)%path_and_outcome_1621_PCL = &
         				processing_information(i,j)%path_and_outcome_1621_PCL + good_retrieval  
		endif

	   if (optical_thickness_16_final_PCL(i,j) > 0. .and. effective_radius_16_final_PCL(i,j) > 0.) & 
         processing_information(i,j)%path_and_outcome_16_PCL = &
         					processing_information(i,j)%path_and_outcome_16_PCL + good_retrieval  
      if (optical_thickness_16_final(i,j) > 0. .and. effective_radius_16_final(i,j) > 0.) &
         processing_information(i,j)%path_and_outcome_16 = &
                        processing_information(i,j)%path_and_outcome_16 + good_retrieval

      if( ( c2_sensor_id == OCI_ID ) .or. ( c2_sensor_id == OCIS_ID ) ) then
        if (optical_thickness_22_final_PCL(i,j) > 0. .and. &
          effective_radius_22_final_PCL(i,j) > 0.) &
          processing_information(i,j)%path_and_outcome_22_PCL = &
          processing_information(i,j)%path_and_outcome_22_PCL + good_retrieval
        if (optical_thickness_22_final(i,j) > 0. .and. &
          effective_radius_22_final(i,j) > 0.) &
          processing_information(i,j)%path_and_outcome_22 = &
          processing_information(i,j)%path_and_outcome_22 + good_retrieval
      else
	     if (optical_thickness_37_final_PCL(i,j) > 0. .and. &
          effective_radius_37_final_PCL(i,j) > 0.) & 
          processing_information(i,j)%path_and_outcome_37_PCL = &
          processing_information(i,j)%path_and_outcome_37_PCL + good_retrieval  
         					
	     if (optical_thickness_37_final(i,j) > 0. .and. &
          effective_radius_37_final(i,j) > 0.) & 
          processing_information(i,j)%path_and_outcome_37 = &
          processing_information(i,j)%path_and_outcome_37 + good_retrieval  
		endif
		   
        
       if (processing_information(i,j)%optical_thickness_GC /= 0 .and.  &
           processing_information(i,j)%band_used_for_optical_thickness ==1 )then
         processing_information(i,j)%rayleigh_correction = 1 
       else
         processing_information(i,j)%rayleigh_correction = 0
       endif 


! initialize the variable first of all. -- GW 4.6.05
       processing_information(i,j)%multi_layer_cloud = 0

       if (cloudsummary(i,j)%cloudmask_determined) then

         if (processing_information(i,j)%path_and_outcome == 1)then ! dec.tree stop
            processing_information(i,j)%multi_layer_cloud = 1

         elseif (processing_information(i,j)%path_and_outcome == 2 .or. &
         		 processing_information(i,j)%path_and_outcome == 10 ) then ! water cloud
            if (cloud_layer_flag(i,j) < 2 .or. ml_test_flag(i,j) == 16) then 
               processing_information(i,j)%multi_layer_cloud = 2 ! SL water
            else 
               processing_information(i,j)%multi_layer_cloud = 3 ! ML water
            endif

         elseif(processing_information(i,j)%path_and_outcome == 3 .or. &
         		processing_information(i,j)%path_and_outcome == 11 ) then !ice cloud
            if (cloud_layer_flag(i,j) < 2 .or. ml_test_flag(i,j) == 16) then 
               processing_information(i,j)%multi_layer_cloud = 4 ! SL ice
            else 
               processing_information(i,j)%multi_layer_cloud = 5 ! ML ice
            endif

         elseif (processing_information(i,j)%path_and_outcome == 4 .or. &
         		 processing_information(i,j)%path_and_outcome == 12 ) then ! unknown cloud
            if (cloud_layer_flag(i,j) < 2 .or. ml_test_flag(i,j) == 16) then 
               processing_information(i,j)%multi_layer_cloud = 6 ! SL unknown
            else 
               processing_information(i,j)%multi_layer_cloud = 7 ! ML unknown
            endif
         endif

       else
          processing_information(i,j)%multi_layer_cloud = 0 
       endif

! set the information for CSR QA -- GW. 4.7.05
       processing_information(i,j)%CSR_flag = 0
       processing_information(i,j)%CSR_flag = CSR_flag_array(i,j)

! set the information for the ML Test QA -- GW. 5.13.09
	   processing_information(i,j)%ml_test_mark = 0
	   processing_information(i,j)%ml_test_mark = ml_test_flag(i,j)

     enddo
   enddo
 end subroutine set_quality_data

 end module set_quality_data_module
