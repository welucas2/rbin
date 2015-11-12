!William Lucas (wel2@st-andrews.ac.uk)
!Skeleton program providing basic framework to run over a sequence of
!sphNG output binaries, performing an operation on each.
!Bugs:

!Note that both this program and rbin must be compiled with the same default
!real precision as the code which produced the binaries. This essentially means
!that a flag for default double precision should always be used. Care should then
!be taken when using any other library that expects to be given data of a certain
!precision (e.g. PGPLOT, which will only accept single precision reals, will
!only be able to use reals declared as REAL*4).

!SAMPLE COMPILE COMMANDS (assuming double precision binaries):
! - gfortran -o <executable name> -fdefault-real-8 rbin.f90 read_sphng_skeleton.f90
! - ifort -o <executable name> -autodouble rbin.f90 read_sphng_skeleton.f90

PROGRAM read_skeleton
  USE rbin
  IMPLICIT NONE
  
  !File read variables
  character(len=4) :: fprefix
  character(len=7) :: fname
  integer :: fstart, fend, fstep
  
  integer :: ii !dummy loop variable
  logical :: ifirst !used run things only the first time through the main loop  
  
  !Firstly get the files we're going to be working on.
  write(*,"(A)",advance="no") "Enter binary file prefix: "
  read (*,"(A4)") fprefix
  DO
    write(*,"(A)",advance="no") "Enter file start and end numbers: "
    read *, fstart, fend
    IF (fstart <= 0 .OR. fend < fstart .OR. fend > 999) THEN
      IF (fstart <= 0) print *, "The starting number must be positive."
      IF (fend < fstart) print *, "The ending number must be >= the starting."
      IF (fend > 999) print *, "The ending number must not be greater than 999."
    ELSE
      EXIT
    END IF
  END DO
  !write(*,"(A)",advance="no") "Enter file number step: "
  !read (*,*,iostat=ioerr) fstep
  fstep = 1 !read each file
  
  !!!!!!!!!!!!!!!
  !If creating an output file to contain data from all binaries, OPEN it here.
  !!!!!!!!!!!!!!!

  contiguous = .TRUE. !Read the entire binary and store data in memory as one block.
                      !Only important when working with MPI runs.
                      
  ifirst = .TRUE.
  write (*,"(A)") "Processing data..."
  
  !****    ****    ****    ****
  !MAIN PROGRAM LOOP OVER FILES
  !****    ****    ****    ****
  write (*,"(A)") "Processing data..."
  ifirst = .TRUE.
  DO ii=fstart,fend,fstep
    !Determine file name, and open it.
    IF (ii < 10) THEN
      write(fname,"(A4,A2,I1)") fprefix, "00", ii
    ELSE IF (ii < 100) THEN
      write(fname,"(A4,A1,I2)") fprefix, "0", ii
    ELSE
      write(fname,"(A4,I3)") fprefix, ii
    END IF
    write (*,*) " "
    write (*,"(A,A7)") "Working on file ",fname
    open(10,file=fname,status='unknown',action='read',form='unformatted')
    
    !!!!!!!!!!!!!!!
    !If creating a new output file for each binary, OPEN it here.
    !!!!!!!!!!!!!!!
  
    !read the data from file
    CALL read_sphng_header
    CALL read_sphng_data
    !and close it.
    close(10)
    
    IF (ifirst) THEN
      !Do anything here which only has to be done once, such as defining derived
      !units (energy, momentum, etc.) or noting the gas particle mass.
    END IF
    
    !!!!!!!!!!!!!!!
    !Do work specific to this binary here.
    !!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!
    !Output data specific to this binary here
    !!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!
    !If creating a new output file for each binary, CLOSE it here.
    !!!!!!!!!!!!!!!
    
    !Clean up rbin by deallocating the arrays created in read_sphng_data
    CALL deallocate_all_arrays

    IF (ifirst) ifirst = .FALSE. !No longer first run through!
    
  END DO
  
  !!!!!!!!!!!!!!!
  !If creating an output file to contain data from all binaries, CLOSE it here.
  !!!!!!!!!!!!!!!
  
END PROGRAM read_skeleton

