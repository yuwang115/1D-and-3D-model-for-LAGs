
FUNCTION ReadSMB(Model,nodenumber,VarIn)  RESULT(VarOut)
  
  USE types
  USE SolverUtils
  
  IMPLICIT NONE
  
  TYPE(Model_t) :: Model
  INTEGER       :: nodenumber
  REAL(kind=dp) :: VarIn
  REAL(kind=dp) :: VarOut

  ! Hard coded file information
  INTEGER,PARAMETER           :: FileLength = 1340, FileUnit = 39
  CHARACTER(len=MAX_NAME_LEN) :: FileName='smb_MGGZ_timestep.asc'
  
!  REAL(KIND=dp)               :: tt, dt  
!  TYPE(variable_t), POINTER   :: Timevar
  INTEGER                     :: ii

!  REAL(KIND=dp),SAVE          :: TimeOffset  
  LOGICAL,SAVE                :: FirstTime = .TRUE.
  REAL(kind=dp),ALLOCATABLE,SAVE :: bmb(:)
  
!  Timevar => VariableGet( Model % Variables,'Time')
!  tt = TimeVar % Values(1)
!  dt = Model % Solver % dt

  IF (FirstTime) THEN
!    TimeOffset = tt - dt
    ALLOCATE(bmb(FileLength))
    OPEN(UNIT=FileUnit, FILE=FileName)
    DO ii = 1,FileLength
      READ(FileUnit,*) bmb(ii)
    END DO
    CLOSE(UNIT=FileUnit)
    FirstTime = .FALSE.
  END IF

!  tt = tt - TimeOffset
!  ii = floor(tt/dt)

  ii = floor(VarIn)
  
  IF (ii.GT.FileLength) CALL FATAL("ReadBMB","File too short...")

  VarOut = bmb(ii)

END FUNCTION ReadSMB
