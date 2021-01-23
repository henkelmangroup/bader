  MODULE automatedcp_mod
    USE kind_mod
    USE matrix_mod
    USE bader_mod
    USE charge_mod 
    USE options_mod
    USE ions_mod
    USE io_mod
    USE ions_mod
    USE weight_mod
    USE dsyevj3_mod
    USE critpoint_mod
    IMPLICIT NONE
    CONTAINS 

    SUBROUTINE AutoCP(bdr,chg,opts,ions,stat)
      TYPE (bader_obj) :: bdr
      TYPE (charge_obj) :: chg
      TYPE (ions_obj)  :: ions
      TYPE (options_obj) :: opts
      INTEGER :: stat
      LOGICAL :: isExhausted
      !stat 1 means pass 0 means fail
      isExhausted = .FALSE.
      DO WHILE ( .NOT. isExhausted .AND. stat /= 1) 
        PRINT *, "Calling critpoint_find"
        CALL critpoint_find(bdr,chg,opts,ions,stat)
        PRINT *, "Stat came back as ", stat
        isExhausted = CheckExhaustion(opts)
        PRINT *, "isExhausted came back as ", isExhausted
        IF (stat /= 1) THEN
          PRINT *, "Adjusting parameters to be stricter"
          CALL AdjustParameters(opts)
        END IF
      END DO
    END SUBROUTINE AutoCP

    FUNCTION CheckExhaustion(opts)
      ! This function checks if the parameter options are exhausted
      TYPE (options_obj) :: opts
      LOGICAL :: CheckExhaustion
      CheckExhaustion = .FALSE.
      PRINT *, "In check exhaustion"
      IF (opts%par_sr<=2 .AND. opts%par_newtonr <= 0.000000001 .AND. &
          opts%par_gradfloor <= 0.000000001) THEN
        PRINT *, "Parameters are strict enough"
        PRINT *, "enableDensityDescend is ", opts%enableDensityDescend
        PRINT *, "enableCHGCARSmoothening is ", opts%enableCHGCARSmoothening
        IF (opts%enableDensityDescend  .AND. &
            opts%enableCHGCARSmoothening ) THEN
          CheckExhaustion = .TRUE.
          PRINT *, "The strictest setting has been used."
          PRINT *, "enableDensityDescend is ", opts%enableDensityDescend
          PRINT *, "enableCHGCARSmoothening is ", opts%enableCHGCARSmoothening
        ELSE IF (.NOT. opts%enableDensityDescend ) THEN
          PRINT *, "Enabling DensityDescend for minima finding"
          PRINT *, "Setting other parameters back to default values"
          opts%enableDensityDescend = .TRUE.
          ! In case parameters has been unnecessarily strict, reset everything
          opts%par_newtonr = 0.0001
          opts%par_gradfloor = 0.0001
          opts%par_sr = 2
        ELSE IF (.NOT. opts%enableCHGCARSmoothening ) THEN
          PRINT *, "Enabling CHGCAR Smoothening"
          PRINT *, "WARNING :: This function is only designed to work with &
            AFLOW"
          PRINT *, "Setting other parameters back to default values"
          opts%enableCHGCARSmoothening = .TRUE.
          opts%par_newtonr = 0.0001
          opts%par_gradfloor = 0.0001
          opts%par_sr = 2
        END IF

      END IF
      RETURN 
    END FUNCTION checkExhaustion

    SUBROUTINE AdjustParameters(opts)
      TYPE (options_obj) :: opts
      IF (opts%par_newtonr > 0.000000001) THEN
        opts%par_newtonr = opts%par_newtonr * 0.1
      END IF
      IF ( opts%par_gradfloor > 0.000000001) THEN
        opts%par_gradfloor = opts%par_gradfloor * 0.1
      END IF
      IF ( opts%par_gradfloor <= 0.000000001 .AND. &
           opts%par_newtonr <= 0.000000001 .AND. opts%par_sr > 2) THEN
        opts%par_sr = opts%par_sr - 1
      END IF
      PRINT *, "Search failed to satisfy constrains! Adjusting parameters:"
      PRINT *, "parnewtonr is now ", opts%par_newtonr
      PRINT *, "pargradfloor is now ", opts%par_gradfloor
      PRINT *, "parsr is now ", opts%par_sr
    END SUBROUTINE AdjustParameters
 
  END MODULE
