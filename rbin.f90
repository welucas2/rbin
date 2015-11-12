!William Lucas, University of St Andrews
!Module contains variables to store data read in from an sphng binary dump.
!Subroutine to read data in is also here.
!This is a transliteration of Ian's read code rather than rdump.
!Created: 30 November 2012
!Last modified: 18 Feburary 2015
!Revision list:
!     - 18 Feb 2015 Enable reading of tagged files
!Bugs: - There seems to be an issue in OSX wherein deallocated arrays are still retained
!        in inactive memory. If reading a lot of files, this can entirely deplete memory.
!Future: Reuse arrays - only reallocate if necessary. Should prevent this leak on OSX.

!IMPORTANT:
!  If compiled with ifort there may be issues with reading large files.
!  If on linux, unlimiting the stacksize first should help; on OSX, where this
!  isn't possible, compile with -heap-arrays 6400 to force creation of temporary
!  arrays on the heap instead. (This is just a number I've found to work :)
!  gfortran is by default more open to using the heap, so this isn't a problem.

!  If working on a multiple-block (MPI) file, choose 'contiguous' before use.
!  With .TRUE. the whole file is read and all data stored in one array.
!  With .FALSE. a single block is read and control passed back to the
!  controlling subroutine - any data already in the arrays is destroyed.

!  Ensure module is scoped correctly or data may be lost - if you don't
!  feel like being careful, use global scope ('USE rbin' in the main program).

MODULE rbin
IMPLICIT NONE

!Values to use with gradh, alpha, MHD, extrernal forces, and planetesimal bits
integer, parameter :: igradh = 1 !0 = off, 1 = on
integer, parameter :: imhd = 0 !0 = off, 1 = on
integer, parameter :: iexf = 0 !external force
integer, parameter :: imigrate = 0 !0 = off, 1 = on
real :: pmrate, rorbitmax

!NTAB is the number of entries in iv which stores the state of the RNG
integer, parameter :: NTAB = 32

logical :: contiguous
logical :: tagged
integer :: iblock !iblock is the block we're currently working on

integer :: nparttot !total number of particles
integer, allocatable :: npartblocks(:) !number of particles in each block

integer*8 :: iuniquemax

!npart contains number of current data, whether for 1 block or a whole simulation
integer :: npart, n1, n2, nreassign, naccrete, nkill, nblocks, iyr, idum, iv(NTAB)
integer :: iplanetesimals, irotpot, idragscheme
integer :: numberarray, icount, icountsink
integer :: nptmass

real :: time, dtmaxdp, gamma, rhozero, RK2, escap, tkin, tgrav, tterm
real :: anglostx, anglosty, anglostz, specang, ptmassin
real :: tmag, Bextx, Bexty, Bextz
real :: hzero, uzero_n2, hmass, gapfac, sdprof
real :: rorbit_orig, min_rplan, max_rplan, planetesimalmass, coremass_orig, coremass


real :: udist, umass, utime, umagfd

integer, allocatable :: isteps(:)
integer*8, allocatable :: iunique(:)
integer*1, allocatable :: iphase(:)
real, allocatable :: xyzmh(:,:), vxyzu(:,:), rho(:)
real*4, allocatable :: dgrav(:), gradh(:), gradhsoft(:), alphaMM(:), poten(:)

!We can't allocate sink arrays at runtime because number of sinks is only given per block;
!contiguous-mode read wouldn't work without moving data between old and new arrays, which
!could be costly.
integer, parameter :: nsinkmax = 2000
integer :: listpm(nsinkmax)
real :: spinx(nsinkmax), spiny(nsinkmax), spinz(nsinkmax)
real :: angaddx(nsinkmax), angaddy(nsinkmax), angaddz(nsinkmax)
real ::  spinadx(nsinkmax), spinady(nsinkmax), spinadz(nsinkmax)

!RT data
integer*2, allocatable :: nneigh(:)
real, allocatable :: e(:), rkappa(:), cv(:), rlambda(:), edd(:)
real, allocatable :: force(:,:)
real*4, allocatable :: dlnTdlnP(:), adiabaticgradient(:)
real*4, allocatable :: pressure(:,:), viscosity(:,:), gravity(:,:), radpres(:,:)

private :: read_sphng_block, allocate_arrays, deallocate_arrays, reallocate_arrays, &
  allocate_arrays_RT, deallocate_arrays_RT, reallocate_arrays_RT

CONTAINS

  !This subroutine reads the data at the top of the binary file and stores useful
  !data as a module global variable.
  SUBROUTINE read_sphng_header
  IMPLICIT NONE
  
  integer :: i !dummy variable
  integer :: ierr !error return from extract
  
  character(len=100) :: fileident
  integer :: number
  
  character*16, dimension(128) :: tagsreal
  
  integer, parameter :: idimhead = 30
  real, dimension(idimhead) :: rheader
  
  print *, "Starting file read..."
  
  read(10)
  read(10) fileident
  print *, trim(fileident)
  IF (fileident(1:1).NE.'F') THEN
    print *, "File does not contain a full dump, exiting!"
    CLOSE(10)
    STOP
  END IF
  IF (fileident(2:2).EQ.'T') THEN
    tagged = .TRUE.
    print *, "File is tagged."
  ELSE
    tagged = .FALSE.
    print *, "File is not tagged."
  ENDIF
  
  read(10) number
  IF (tagged) read(10) !skip tags here
  IF (number == 6) THEN
    read(10) nparttot,n1,n2,nreassign,naccrete,nkill
    nblocks = 1
  ELSE IF (number == NTAB + 12) THEN
    read(10) nparttot,n1,n2,nreassign,naccrete,nkill,nblocks,iyr,&
      idum,(iv(i),i=1,NTAB),iplanetesimals,irotpot,idragscheme
  ELSE
    read(10) nparttot,n1,n2,nreassign,naccrete,nkill,nblocks
  END IF
  write (*,*) 'nparttot = ',nparttot,' nblocks ',nblocks
  allocate( npartblocks(nblocks) )
  
  DO i = 1, 3
    read(10) 
  END DO
  read(10) number
  IF (number == 1) THEN
    IF (tagged) read(10)
    read(10) iuniquemax
  ELSE
    iuniquemax = nparttot
  END IF
  print *, "iuniquemax =", iuniquemax
  
  read(10) number
  print *, "number real", number
  
  IF (tagged) THEN
    write(*,"(A,I3)") "Reading tags from file, number =", number
    read(10) tagsreal(1:number)
  ELSE
    print *, "Simulating tags for non-tagged file."
    tagsreal(1:30) = &
      (/'gt              ','dtmax           ','gamma           ', &
      'rhozero         ','RK2             ','escap           ', &
      'tkin            ','tgrav           ','tterm           ', &
      'anglostx        ','anglosty        ','anglostz        ', &
      'specang         ','ptmassin        ','tmag            ', &
      'Bextx           ','Bexty           ','Bextz           ', &
      'hzero           ','uzero_n2        ','hmass           ', &
      'gapfac          ','                ','sdprof          ', &
      'rorbit_orig     ','min_rplan       ','max_rplan       ', &
      'planetesimalmass','coremass_orig   ','coremass        '/)
  END IF  
  
  read(10) (rheader(i),i=1,min(number,idimhead))
  
  CALL extract('gt',time,rheader,tagsreal,number,ierr)
  IF (ierr.NE.0) STOP
  print *, "Time is", time
  CALL extract('dtmax',dtmaxdp,rheader,tagsreal,number,ierr)
  IF (ierr.NE.0) STOP
  CALL extract('gamma',gamma,rheader,tagsreal,number,ierr)
  IF (ierr.NE.0) STOP
  CALL extract('rhozero',rhozero,rheader,tagsreal,number,ierr)
  IF (ierr.NE.0) STOP
  CALL extract('RK2',RK2,rheader,tagsreal,number,ierr)
  IF (ierr.NE.0) STOP
  CALL extract('escap',escap,rheader,tagsreal,number,ierr)
  CALL extract('tkin',tkin,rheader,tagsreal,number,ierr)
  CALL extract('tgrav',tgrav,rheader,tagsreal,number,ierr)
  CALL extract('tterm',tterm,rheader,tagsreal,number,ierr)
  CALL extract('anglostx',anglostx,rheader,tagsreal,number,ierr)
  CALL extract('anglosty',anglosty,rheader,tagsreal,number,ierr)
  CALL extract('anglostz',anglostz,rheader,tagsreal,number,ierr)
  CALL extract('specang',specang,rheader,tagsreal,number,ierr)
  CALL extract('ptmassin',ptmassin,rheader,tagsreal,number,ierr)
  IF (imhd == 1) THEN
    CALL extract('tmag',tmag,rheader,tagsreal,number,ierr)
    CALL extract('Bextx',Bextx,rheader,tagsreal,number,ierr)
    CALL extract('Bexty',Bexty,rheader,tagsreal,number,ierr)
    CALL extract('Bextz',Bextz,rheader,tagsreal,number,ierr)
    print *, 'External field found, Bext = ',Bextx,Bexty,Bextz
  END IF
  CALL extract('hzero',hzero,rheader,tagsreal,number,ierr)
  IF (iexf == 10 .AND. ierr /= 0) STOP      
  CALL extract('uzero_n2',uzero_n2,rheader,tagsreal,number,ierr)
  IF (n2 > 0) THEN
    print *, ' read u for surrounding medium = ',uzero_n2
  ENDIF
  CALL extract('hmass',hmass,rheader,tagsreal,number,ierr)
  CALL extract('gapfac',gapfac,rheader,tagsreal,number,ierr)
  CALL extract('sdprof',sdprof,rheader,tagsreal,number,ierr)
  IF (ierr /= 0) THEN
    sdprof = -0.5
    print *, 'Surface density pre vary (goes as r^-0.5)'
  ENDIF
  CALL extract('rorbit_orig',rorbit_orig,rheader,tagsreal,number,ierr)
  IF (imigrate > 0) rorbitmax = (rorbitmax - rorbit_orig)/pmrate
  CALL extract('min_rplan',min_rplan,rheader,tagsreal,number,ierr)
  CALL extract('max_rplan',max_rplan,rheader,tagsreal,number,ierr)
  IF (ierr /= 0) max_rplan = 1.0E6
  print *, 'Setting planetesimal boundaries; check values.'
  print *, 'r_min = ', min_rplan
  print *, 'r_max = ', max_rplan
  CALL extract('planetesimalmass',planetesimalmass,rheader,tagsreal,number,ierr)
  print *, 'Planetesimal mass ', planetesimalmass
  CALL extract('coremass_orig',coremass_orig,rheader,tagsreal,number,ierr)
  print *, 'Core mass orig ', coremass_orig
  CALL extract('coremass',coremass,rheader,tagsreal,number,ierr)
  print *, 'Core mass running record ', coremass

  read(10) number
  read(10) number
  IF (number < 3) THEN
    print *, "Error in rbin, nreal8 too small in header section"
    STOP
  END IF
  IF (tagged) read(10) !skip tags
  IF (imhd == 1) THEN
    IF (number > 3) THEN
      read(10) udist, umass, utime, umagfd
    ELSE
      print *, 'WARNING: no mag field units in rdump'
      read(10) udist, umass, utime
      !umagfdi = umagfd
      !Will need to define B field units another way      
    END IF
  ELSE
    read(10) udist, umass, utime
  END IF

!Old code to read general dump info
!  read(10) time, dtmaxdp, gamma, rhozero, RK2, escap, tkin, tgrav, tterm, &
!    anglostx, anglosty, anglostz, specang, ptmassin
!  print *, "Time is", time
!  read(10)
!  read(10)
  
!  read(10) udist, umass, utime
  print *, "Distance, mass, time units are:",udist,umass,utime,"in cgs"
  
  read(10) number
  numberarray = number/nblocks
  print *, "Array types", number, numberarray
  !IF (numberarray /= 2) THEN
  !  print *, "Expected 2 for no MHD & no RT - exiting"
  !  CLOSE(10)
  !  STOP
  !END IF
  
  icount = 0
  icountsink = 0
  
  !At this point we are onto the first data block
  iblock = 1
  
  RETURN
  END SUBROUTINE read_sphng_header
  
  !Reads data from the binary. If binary is from an MPI-enabled simulation (ie. is
  !separated into blocks), it will call read_sphng_block multiple times to form a
  !contiguous set of data if the module variable 'contiguous' is set .TRUE. If it is
  !.FALSE., read_sphng_block will be called only once, and any data already in the
  !arrays will be replaced with new entries.
  SUBROUTINE read_sphng_data
  IMPLICIT NONE
  integer :: ii
  !If only 1 block exists or non-contiguous mode, only read once.
  IF (.NOT. contiguous .OR. nblocks == 1) THEN
    CALL read_sphng_block
  !If we want an entire MPI run, call multiple times.
  ELSE IF (contiguous) THEN
    DO ii = 1, nblocks
      CALL read_sphng_block
    END DO
  ELSE
    print *, "Unsure how to read data:"
    print *, "nblocks = ", nblocks
    IF (contiguous) THEN
      print *, "contiguous mode: ON"
    ELSE
      print *, "contiguous mode: OFF"
    END IF
    print *, "Exiting..."
    STOP
  END IF
  RETURN
  END SUBROUTINE read_sphng_data
  
  !Reads in a single block of data.
  SUBROUTINE read_sphng_block
  IMPLICIT NONE
  
  character(len=16) :: tagi
  
  integer :: i !dummy integer
  integer*8 :: number8
  integer, dimension(8) :: nums, numssink, numsRT
  integer :: ic1, ic2, ic3
  
  print *, " "
  print *, " "
  print *, "Reading block ", iblock
  
  read(10) number8, nums(1:8)
  print *, nparttot, nblocks, iblock, number8, number8*nblocks
  npart = number8
  npartblocks(iblock) = npart
  
  print *, "Block contains ", npart, "particles, current icount ", icount
  print *, "nums ", nums(1:8)
  
  read(10) number8, numssink(1:8)
  nptmass = number8
  print *, "numssink ", numssink(1:8)
  IF (numberarray == 3) THEN
    read(10) number8, numsRT(1:8)
    print *, "numsRT ", numsRT(1:8)
  END IF
  
  ic1 = icount + 1
  ic2 = icount + 2
  ic3 = icount + 3
  
  !Re/allocate as needed - isteps is taken to be representative of all arrays
  IF (.NOT. allocated(isteps)) THEN
    CALL allocate_arrays
  !Need to reallocate if we're on blockwise read or if we're starting to read a new file
  ELSE IF (.NOT. contiguous .OR. iblock == 1) THEN
    CALL reallocate_arrays
  END IF
  
  IF (tagged) read(10) !skip tag
  read(10) isteps(icount+1:icount+npart)
  print *, "isteps: ", isteps(ic1:ic3)
  IF (nums(1) >= 2) THEN
    IF (tagged) read(10) !skip tag
    read(10) !skip reading listinactive
  END IF
  IF (tagged) THEN
    read(10) tagi
    IF (trim(tagi) /= "iphase") THEN
      print*," WARNING: iphase does not match tag "//trim(tagi)
    END IF
  END IF
  read(10) iphase(icount+1:icount+npart)
  print *, "iphase: ", iphase(ic1:ic3)
  
  !In rdump, only if nums1(5) >= 1
  IF (nums(5) >= 1) THEN
    IF (tagged) read(10)
    read(10) iunique(icount+1:icount+npart)
  END IF
  
  print *, "Reading xyzmh..."
  DO i = 1, 5
    IF (tagged) read(10) tagi !skip tag
    read(10) xyzmh(i,icount+1:icount+npart)
  END DO
  print *, "x: ", xyzmh(1,ic1:ic3)
  print *, "m: ", xyzmh(4,ic1:ic3)
  print *, "h: ", xyzmh(5,ic1:ic3)
  
  print *, "Reading vxyzu..."
  DO i = 1, 4
    IF (tagged) read(10) tagi !skip tag
    read(10) vxyzu(i,icount+1:icount+npart)
  END DO
  print *, "vx: ", vxyzu(1,ic1:ic3)
  print *, "u: ", vxyzu(4,ic1:ic3)
  
  print *, "Reading rho..."
  IF (tagged) read(10) tagi !skip tag
  read(10) alphaMM(icount+1:icount+npart)
  rho(icount+1:icount+npart) = alphaMM(icount+1:icount+npart)
  print *, "rho: ", rho(ic1:ic3)
  
  print *, "nums(7) is ", nums(7)
  IF (igradh == 1) THEN
    IF (nums(7) >= 2) THEN
      print *, "Reading gradh..."
      IF (tagged) read(10) tagi !skip tag
      print *, "tag is ", trim(tagi)  
      read(10) gradh(icount+1:icount+npart)
      print *, "gradh: ", gradh(ic1:ic3)
    END IF
    IF (nums(7) >= 3) THEN
      print *, "Reading gradhsoft..."
      IF (tagged) read(10) tagi !skip tag
      print *, "tag is ", trim(tagi)  
      read(10) gradhsoft(icount+1:icount+npart)
      print *, "gradhsoft: ", gradhsoft(ic1:ic3)
    END IF
  ELSE
    print *, "Reading dgrav..."
    IF (tagged) read(10) tagi !skip tag
    !print *, "tag is ", trim(tagi)
    read(10) dgrav(icount+1:icount+npart)
    print *, "dgrav: ", dgrav(ic1:ic3)  
  END IF 
  
  
!  print *, "nums(7) is ", nums(7)
!  IF (tagged) read(10) tagi !skip tag
!  DO i = 1, nums(7)-2
!    read(10)
!  END DO

  !Old non-gradh dumps may not include alphaMM
  IF (nums(7) >= 4 .OR. igradh == 0 .AND. nums(7) >= 3) THEN 
    print *, "Reading alphaMM..."
    IF (tagged) read(10) tagi !skip tag
    read(10) alphaMM(icount+1:icount+npart)
    print *, "alphaMM: ", alphaMM(ic1:ic3)
  END IF

  IF (nums(7) >= 5) THEN
    IF (tagged) read(10) tagi !skip tag
    read(10) poten(icount+1:icount+npart)
    print *, "poten: ", poten(ic1:ic3)  
  END IF
  
  !If nptmass is zero, this may be incorrect (though it runs anyway).
  !It may be safer to read to a buffer and only store if nptmass > 0.
  print *, "Reading sink particle data, nptmass: ", nptmass
  IF (tagged) read(10) tagi !skip tag
  read(10) listpm(icountsink+1:icountsink+nptmass)
  IF (tagged) read(10) tagi !skip tag
  read(10) spinx(icountsink+1:icountsink+nptmass)
  IF (tagged) read(10) tagi !skip tag
  read(10) spiny(icountsink+1:icountsink+nptmass)
  IF (tagged) read(10) tagi !skip tag
  read(10) spinz(icountsink+1:icountsink+nptmass)
  IF (tagged) read(10) tagi !skip tag
  read(10) angaddx(icountsink+1:icountsink+nptmass)
  IF (tagged) read(10) tagi !skip tag
  read(10) angaddy(icountsink+1:icountsink+nptmass)
  IF (tagged) read(10) tagi !skip tag
  read(10) angaddz(icountsink+1:icountsink+nptmass)
  IF (tagged) read(10) tagi !skip tag
  read(10) spinadx(icountsink+1:icountsink+nptmass)  
  IF (tagged) read(10) tagi !skip tag
  read(10) spinady(icountsink+1:icountsink+nptmass)  
  IF (tagged) read(10) tagi !skip tag
  read(10) spinadz(icountsink+1:icountsink+nptmass)
  
!I'm just not sure about this point onwards - if RT comes into things
! (and it probably won't) then check it. But that won't be the case.
  
  print *, "numssink(6) = ", numssink(6)
  DO i = 1, numssink(6)-9
    read(10)
  END DO
  DO i = 1, numssink(8)
    IF (tagged) read(10) !skip tag
    read(10)
  END DO
  
  !Read RT data if it's there
  !NB Original indices were from 1 to npart; I'm guessing they should be changed.
  !Also NB - is this old? Some parts don't work logically...
  IF (numberarray == 3) THEN
    print *, "Reading RT", numberarray, numsRT(1:8)
    IF (.NOT.allocated(nneigh)) THEN
      CALL allocate_arrays_RT
    ELSE IF (.NOT. contiguous .OR. iblock == 1) THEN
      CALL reallocate_arrays_RT
    END IF
    IF (numsRT(3) == 1) read(10) nneigh(icount+1:icount+npart)
    read(10) e(icount+1:icount+npart)
    read(10) rkappa(icount+1:icount+npart)
    read(10) cv(icount+1:icount+npart)
    read(10) rlambda(icount+1:icount+npart)
    read(10) edd(icount+1:icount+npart)
    IF (numsRT(6) == 8) THEN
      DO i = 1, 3
        read(10) force(i,icount+1:icount+npart)
      END DO
    END IF
    IF (numsRT(7) == 1) read(10) dlnTdlnP(icount+1:icount+npart)
    IF (numsRT(7) == 13 .OR. numsRT(7) == 14) THEN
      read(10) dlnTdlnP(icount+1:icount+npart)
      IF (numsRT(7) == 11) read(10) adiabaticgradient(icount+1:icount+npart)
      DO i = 1, 3
        read(10) pressure(i,icount+1:icount+npart)
      END DO
      DO i = 1, 3
        read(10) viscosity(i,icount+1:icount+npart)
      END DO
      DO i = 1, 3
        read(10) gravity(i,icount+1:icount+npart)
      END DO
      DO i = 1, 3
        read(10) radpres(i,icount+1:icount+npart)
      END DO
    END IF
  END IF
  
  IF (contiguous) THEN
    icount = icount + npart
    icountsink = icountsink + nptmass
  ELSE
    icount = 0
    icountsink = 0
  END IF
  iblock = iblock + 1
  
  RETURN
  END SUBROUTINE read_sphng_block
  
!---------------------------------------------------------------
!extract subroutine taken directly from rdump and reformatted to fit in
!with the rest of rbin.
  SUBROUTINE extract(tag,rval,rarr,tags,ntags,ierr)
  character(len=*), intent(in)  :: tag
  real, intent(out) :: rval
  real, intent(in)  :: rarr(:)
  character(len=16), intent(in)  :: tags(:)
  integer, intent(in)  :: ntags
  integer, intent(out) :: ierr
  logical :: matched
  integer :: i

  ierr = 1
  matched = .FALSE.
  rval = 0.0 ! default if not found
  over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
      if (size(rarr) >= i) then
        rval = rarr(i)
          matched = .true.
      endif
      exit over_tags  ! only match first occurrence
    endif
  enddo over_tags
  if (matched) ierr = 0
  if (ierr.NE.0) print "(a)", &
    ' ERROR: could not find '//trim(adjustl(tag))//' in header'
  END SUBROUTINE extract
  
  
!---------------------------------------------------------------
!These create and destroy all the module's arrays, except the sink arrays and npartblocks.
  
  SUBROUTINE allocate_arrays
  IMPLICIT NONE
  integer :: nalloc
  IF (contiguous .AND. nblocks > 1) THEN
    nalloc = nparttot
  ELSE
    nalloc = npart
  END IF
  print *, "ALLOCATING HYDRO ARRAYS"
  allocate( isteps(nalloc), iphase(nalloc), iunique(nalloc) )
  allocate( xyzmh(5,nalloc), vxyzu(4,nalloc), alphaMM(nalloc), poten(nalloc), rho(nalloc) )
  IF (igradh == 1) THEN
    allocate( gradh(nalloc), gradhsoft(nalloc) )
  ELSE
    allocate( dgrav(nalloc) )
  END IF
  RETURN
  END SUBROUTINE allocate_arrays
  
  SUBROUTINE deallocate_arrays
  IMPLICIT NONE
  print *, "DEALLOCATING HYDRO ARRAYS"
  deallocate( isteps, iphase, iunique )
  deallocate( xyzmh, vxyzu, alphaMM, poten, rho )
  IF (igradh == 1) THEN
    deallocate( gradh, gradhsoft )
  ELSE
    deallocate( dgrav )
  END IF
  RETURN
  END SUBROUTINE deallocate_arrays
  
  SUBROUTINE reallocate_arrays
  IMPLICIT NONE
  CALL deallocate_arrays
  CALL allocate_arrays
  RETURN
  END SUBROUTINE reallocate_arrays
  
  SUBROUTINE allocate_arrays_RT
  IMPLICIT NONE
  integer :: nalloc
  IF (contiguous .AND. nblocks > 1) THEN
    nalloc = nparttot
  ELSE
    nalloc = npart
  END IF
  print *, "Allocating RT arrays"
  allocate( nneigh(nalloc) )
  allocate( e(nalloc), rkappa(nalloc), cv(nalloc), rlambda(nalloc), edd(nalloc) )
  allocate( force(3,nalloc) )
  allocate( dlnTdlnP(nalloc), adiabaticgradient(nalloc) )
  allocate( pressure(3,nalloc), viscosity(3,nalloc), gravity(3,nalloc), radpres(3,nalloc) )
  RETURN
  END SUBROUTINE allocate_arrays_RT
  
  SUBROUTINE deallocate_arrays_RT
  IMPLICIT NONE
  print *, "Deallocating RT arrays"
  deallocate( nneigh )
  deallocate( e, rkappa, cv, rlambda, edd )
  deallocate( force )
  deallocate( dlnTdlnP, adiabaticgradient )
  deallocate( pressure, viscosity, gravity, radpres )
  RETURN
  END SUBROUTINE deallocate_arrays_RT
  
  SUBROUTINE reallocate_arrays_RT
  IMPLICIT NONE
  CALL deallocate_arrays_RT
  CALL allocate_arrays_RT
  RETURN
  END SUBROUTINE reallocate_arrays_RT
  
  !Deallocate all arrays held in this module, including npartblocks.
  SUBROUTINE deallocate_all_arrays
  IMPLICIT NONE
  print *, "Deallocating all arrays:"
  IF (allocated(isteps)) CALL deallocate_arrays
  IF (allocated(nneigh)) CALL deallocate_arrays_RT
  IF (allocated(npartblocks)) deallocate(npartblocks)
  RETURN
  END SUBROUTINE deallocate_all_arrays

END MODULE rbin
