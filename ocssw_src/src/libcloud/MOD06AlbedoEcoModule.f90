module MOD06AlbedoEcoModule

! HI
!****************************************************************************
  ! !F90
  !
  ! !Description:
  !    This module contains the subroutines used to provide albedo and ecosystem
  !     maps for processing a MOD06 granule.
  !    There is only one callable routine, getAlbedoEco, which returns albedo 
  !     values and IGBP ecosystem classification for each pixel of the MOD06 granule,
  !     and for the wavelengths specified.  Also returned are values of snow albedos
  !     by ecosystem class for the wavelengths specified.
  !  
  ! !Callable routines:
  !    getAlbedoEco()
  ! 
  ! !Revision History:
  !
  ! Revision 1.0  2003/12/18  12:43:43  EGMoody
  ! Initial revision.
  !
  ! !Team-Unique Header:
  !   Cloud Retrieval Group, NASA Goddard Space Flight Center
  !
  ! !References and Credits:
  !   Written by
  !    Eric Moody
  !    Climate and Radiation Branch, Code 913
  !    NASA/GSFC
  !    Greenbelt MD 20771
  !    Eric.Moody@gsfc.nasa.gov
  !
  ! !Design Notes:
  !
  ! !END
  !*****************************************************************************

  !Dependencies:

  use nonscience_parameters
  implicit none

  include 'netcdf.inc'

  private

  ! WDR public :: GetNISEType, ReadSnowAlbStats, init_NISE_processing
  public :: ReadSnowAlbStats



  !local variables:
  !counters:
  integer                          :: i,j,k,l,m,n

  !HDF variables:
  integer                           :: HDFstatus
  integer, dimension(10)                         :: hdfStart, hdfStride, hdfEdge
  integer                          :: sds_id, sds_index
  !  HDFstatus: Used for error checking.
  !  hdfstart : Follows HDF conventions. The array element from which to begin reading the
  !              HDF array. Zero based.
  !  hdfstride: Follows HDF conventions. The frequency with which data should be read in
  !              each dimension. Stride 1 means read every element.
  !  hdfedge  : Follows HDF conventions. The number of elements to read in each dimension. 
  !              If start(:) == 0 and stride(:) == 1, setting edge equal to the shape of 
  !              the underlying array reads all elements. 


  INTEGER, PARAMETER                              :: gridsize = 721  
  integer*2, DIMENSION(gridsize, gridsize) :: niseNorth, niseSouth


contains

!	subroutine init_NISE_processing(nise_file)
!	
!		character(len=*), intent(in) :: nise_file
!		integer :: iret
!	
!		iret = Read_NISE(nise_file, gridsize, niseNorth, niseSouth)
!
!
!	
!	end subroutine init_NISE_processing



  ! ----------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------

  subroutine ReadSnowAlbStats  ( StatsFN, NumSnowTypes, NumAlbBnds, numEco, &
                                 AlbedoMean,  errorLevel                    )
                                 
    character (len = *),                 intent( in)  :: StatsFN
    real      ,  &
              dimension(:,0:,:),         intent(inout):: AlbedoMean
    integer   , intent( in)  :: NumAlbBnds,       &
                                                         numEco,           &
                                                         NumSnowTypes
    INTEGER,            INTENT(OUT) :: errorLevel
                                                         
    ! !Description:
    !   This routine will output only a single band of albedo, for the specified amount
    !    of data, to the specified position within the hdf file.
    !  
    ! !Input Parameters:
    !    StatsFN : The name, only, of the new HDF file.
    !    Albedo          : Contains the single band of albedo to be stored.
    !    NumSnowTypes    : 2, Dry or wet, via NISE classes.
    !    NumAlbBnds      : Number of Albedo wavelengths.
    !    numEco          : Number of Ecosystem Classes.
    !
    ! !Output Parameters:
    !    AlbedoMean      : Mean statistic
    !
    ! !Revision History:
    !   See Module revision history at the beginning of the file.
    !
    ! !Team-Unique Header:
    !   Cloud Retrieval Group, NASA Goddard Space Flight Center
    !
    ! !References and Credits:
    !   Written by
    !    Eric Moody
    !    Climate and Radiation Branch, Code 913
    !    NASA/GSFC
    !    Greenbelt MD 20771
    !    Eric.Moody@gsfc.nasa.gov
    !
    ! !Design Notes:
    !
    ! !END
    
    !local variables: 
    !HDF variables:
    integer               :: status
    integer , dimension(10) :: hdfStart, hdfStride, hdfEdge
    integer               :: sds_id, newHDFID, sds_index
    character(len =  200)                :: SDSName
    real      ,  &
              dimension(1:2,1:NumAlbBnds,0:numEco,1:1,1:2)   :: DummyAlb
    !  status    : Used for error checking.
    !  hdfstart  : Follows HDF conventions. The array element from which to begin reading the
    !               HDF array. Zero based.
    !  hdfstride : Follows HDF conventions. The frequency with which data should be read in
    !               each dimension. Stride 1 means read every element.
    !  hdfedge   : Follows HDF conventions. The number of elements to read in each dimension. 
    !               If start(:) == 0 and stride(:) == 1, setting edge equal to the shape of 
    !               the underlying array reads all elements.
    !  sds_id,
    !   newHDFID : HDF SDS ID.
    !  sds_index : Index of the SDS in the HDF file.
    ! SDSName    : Stores the name of the SDS being procesed.

   integer :: FAIL = -1

  !Set error level
  errorLevel = success 

    !************************************************************************************
    ! Open the input file:
    !************************************************************************************
    !Open the albedo HDF file:
    status = nf_open( trim(StatsFN), NF_NOWRITE, newHDFID )
    call cld_fchk( status, __FILE__, __LINE__ )

    !  WDR the errorLevel was passed back and not dealt with.  The above will
    !  halt the program here

    !************************************************************************************
    ! Input the Albedo mean statistic:
    !************************************************************************************   
    !Set up the output hdf variables, in this case we are writing a Belt:
    ! hdfStart, note that we are storing the entire Longitude, Albedo Bands and Ecosystems
    ! so these start at 0, however we are storing a box of Latitude, so this starts at the
    ! starty counter.
    hdfStart ( : ) = 1
    hdfStart ( 5 ) = 1  !Latitude, 0 = NH, 1 = SH
    hdfStart ( 4 ) = 1  !Longitude, only 1, so 0
    hdfStart ( 3 ) = 1  !Eco
    hdfStart ( 2 ) = 1  !AlbBands
    hdfStart ( 1 ) = 1  !NumSnowTypes

    !Strides = 1, since we are not skiping points
    hdfStride( : ) = 1
    
    !Edge is the total number of points being read:
    ! It is the number of Albedo bands, the number of Ecosystems, the Number of Longitude
    ! boxes, and 1 lat box:
    hdfEdge  ( : ) = 1
    hdfEdge  ( 5 ) = 2
    hdfEdge  ( 4 ) = 1
    hdfEdge  ( 3 ) = numEco
    hdfEdge  ( 2 ) = NumAlbBnds
    hdfEdge  ( 1 ) = NumSnowTypes
    
    !Read the Mean data:
    SDSName = 'Snow_Albedo_Year_Mean'
    !determine the sds_index for the SDS:
    status = nf_inq_varid( newHDFID, trim(SDSName), sds_id )
    call cld_fchk( status, __FILE__, __LINE__ )

    !get access to this SDS:
        
    !Read the data:
    status = nf_get_vara_real( newHDFID, sds_id, hdfStart, hdfEdge, DummyAlb )
    call cld_fchk( status, __FILE__, __LINE__ )

    !Store the dummy data in the final array:
    do i = 1, NumAlbBnds
      do j = 0, numEco-1
        do k = 1, NumSnowTypes
           AlbedoMean(i,j,k) = DummyAlb(k,i,j,1,1)
        end do
      end do
    end do

    !Error Checking:
    if (sds_index == fail .or. &
        sds_id    == fail .or. &
        status    == fail)     then
        errorLevel = fail
        return
     end if
    
    !************************************************************************************
    ! Close the file:
    !************************************************************************************
    !Close the  HDF file:
    status = nf_close( newHDFID )
    call cld_fchk( status, __FILE__, __LINE__ )

  end subroutine ReadSnowAlbStats

end module MOD06AlbedoEcoModule
