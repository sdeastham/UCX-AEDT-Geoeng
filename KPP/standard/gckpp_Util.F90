! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Auxiliary Routines File
! 
! Generated by KPP-2.2.1_rs5 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : gckpp_Util.f90
! Time                 : Mon Jun  2 12:32:46 2014
! Working directory    : /home/mpayer/KPP
! Equation file        : gckpp.kpp
! Output root filename : gckpp
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE gckpp_Util

  USE gckpp_Parameters
  IMPLICIT NONE

CONTAINS



! User INLINED Utility Functions

! End INLINED Utility Functions

! Utility Functions from KPP_HOME/util/util
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! UTIL - Utility functions
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ****************************************************************
!                            
! InitSaveData - Opens the data file for writing
!   Parameters :                                                  
!
! ****************************************************************

      SUBROUTINE InitSaveData ()

      USE gckpp_Parameters

      open(10, file='gckpp.dat')

      END SUBROUTINE InitSaveData

! End of InitSaveData function
! ****************************************************************

! ****************************************************************
!                            
! SaveData - Write LOOKAT species in the data file 
!   Parameters :                                                  
!
! ****************************************************************

      SUBROUTINE SaveData ()

      USE gckpp_Global
      USE gckpp_Monitor

      INTEGER i

      WRITE(10,999) (TIME-TSTART)/3600.D0,  &
                   (C(LOOKAT(i))/CFACTOR, i=1,NLOOKAT)
999   FORMAT(E24.16,100(1X,E24.16))

      END SUBROUTINE SaveData

! End of SaveData function
! ****************************************************************

! ****************************************************************
!                            
! CloseSaveData - Close the data file 
!   Parameters :                                                  
!
! ****************************************************************

      SUBROUTINE CloseSaveData ()

      USE gckpp_Parameters

      CLOSE(10)

      END SUBROUTINE CloseSaveData

! End of CloseSaveData function
! ****************************************************************

! ****************************************************************
!                            
! GenerateMatlab - Generates MATLAB file to load the data file 
!   Parameters : 
!                It will have a character string to prefix each 
!                species name with.                                                 
!
! ****************************************************************

      SUBROUTINE GenerateMatlab ( PREFIX )

      USE gckpp_Parameters
      USE gckpp_Global
      USE gckpp_Monitor

      
      CHARACTER(LEN=8) PREFIX 
      INTEGER i

      open(20, file='gckpp.m')
      write(20,*) 'load gckpp.dat;'
      write(20,990) PREFIX
990   FORMAT(A1,'c = gckpp;')
      write(20,*) 'clear gckpp;'
      write(20,991) PREFIX, PREFIX
991   FORMAT(A1,'t=',A1,'c(:,1);')
      write(20,992) PREFIX
992   FORMAT(A1,'c(:,1)=[];')

      do i=1,NLOOKAT
        write(20,993) PREFIX, SPC_NAMES(LOOKAT(i)), PREFIX, i
993     FORMAT(A1,A6,' = ',A1,'c(:,',I2,');')
      end do
      
      CLOSE(20)

      END SUBROUTINE GenerateMatlab

! End of GenerateMatlab function
! ****************************************************************


! End Utility Functions from KPP_HOME/util/util
! End of UTIL function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Shuffle_user2kpp - function to copy concentrations from USER to KPP
!   Arguments :
!      V_USER    - Concentration of variable species in USER's order
!      V         - Concentrations of variable species (local)
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Shuffle_user2kpp ( V_USER, V )

! V_USER - Concentration of variable species in USER's order
  REAL(kind=dp) :: V_USER(NVAR)
! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)

  V(94) = V_USER(1)
  V(88) = V_USER(2)
  V(123) = V_USER(3)
  V(46) = V_USER(4)
  V(114) = V_USER(5)
  V(106) = V_USER(6)
  V(68) = V_USER(7)
  V(55) = V_USER(8)
  V(116) = V_USER(9)
  V(95) = V_USER(10)
  V(1) = V_USER(11)
  V(2) = V_USER(12)
  V(3) = V_USER(13)
  V(4) = V_USER(14)
  V(5) = V_USER(15)
  V(6) = V_USER(16)
  V(7) = V_USER(17)
  V(8) = V_USER(18)
  V(9) = V_USER(19)
  V(10) = V_USER(20)
  V(11) = V_USER(21)
  V(12) = V_USER(22)
  V(13) = V_USER(23)
  V(14) = V_USER(24)
  V(15) = V_USER(25)
  V(16) = V_USER(26)
  V(111) = V_USER(27)
  V(56) = V_USER(28)
  V(96) = V_USER(29)
  V(102) = V_USER(30)
  V(47) = V_USER(31)
  V(120) = V_USER(32)
  V(49) = V_USER(33)
  V(93) = V_USER(34)
  V(57) = V_USER(35)
  V(135) = V_USER(36)
  V(73) = V_USER(37)
  V(122) = V_USER(38)
  V(76) = V_USER(39)
  V(115) = V_USER(40)
  V(104) = V_USER(41)
  V(92) = V_USER(42)
  V(89) = V_USER(43)
  V(77) = V_USER(44)
  V(79) = V_USER(45)
  V(110) = V_USER(46)
  V(124) = V_USER(47)
  V(105) = V_USER(48)
  V(121) = V_USER(49)
  V(61) = V_USER(50)
  V(99) = V_USER(51)
  V(48) = V_USER(52)
  V(129) = V_USER(53)
  V(131) = V_USER(54)
  V(109) = V_USER(55)
  V(130) = V_USER(56)
  V(54) = V_USER(57)
  V(100) = V_USER(58)
  V(63) = V_USER(59)
  V(117) = V_USER(60)
  V(43) = V_USER(61)
  V(136) = V_USER(62)
  V(139) = V_USER(63)
  V(134) = V_USER(64)
  V(137) = V_USER(65)
  V(132) = V_USER(66)
  V(45) = V_USER(67)
  V(81) = V_USER(68)
  V(108) = V_USER(69)
  V(67) = V_USER(70)
  V(39) = V_USER(71)
  V(107) = V_USER(72)
  V(80) = V_USER(73)
  V(64) = V_USER(74)
  V(113) = V_USER(75)
  V(133) = V_USER(76)
  V(127) = V_USER(77)
  V(69) = V_USER(78)
  V(51) = V_USER(79)
  V(52) = V_USER(80)
  V(128) = V_USER(81)
  V(112) = V_USER(82)
  V(119) = V_USER(83)
  V(72) = V_USER(84)
  V(85) = V_USER(85)
  V(62) = V_USER(86)
  V(103) = V_USER(87)
  V(71) = V_USER(88)
  V(53) = V_USER(89)
  V(44) = V_USER(90)
  V(17) = V_USER(91)
  V(18) = V_USER(92)
  V(19) = V_USER(93)
  V(59) = V_USER(94)
  V(138) = V_USER(95)
  V(90) = V_USER(96)
  V(58) = V_USER(97)
  V(60) = V_USER(98)
  V(40) = V_USER(99)
  V(65) = V_USER(100)
  V(36) = V_USER(101)
  V(37) = V_USER(102)
  V(38) = V_USER(103)
  V(20) = V_USER(104)
  V(21) = V_USER(105)
  V(22) = V_USER(106)
  V(23) = V_USER(107)
  V(50) = V_USER(108)
  V(98) = V_USER(109)
  V(78) = V_USER(110)
  V(101) = V_USER(111)
  V(86) = V_USER(112)
  V(118) = V_USER(113)
  V(74) = V_USER(114)
  V(75) = V_USER(115)
  V(84) = V_USER(116)
  V(125) = V_USER(117)
  V(126) = V_USER(118)
  V(91) = V_USER(119)
  V(66) = V_USER(120)
  V(97) = V_USER(121)
  V(87) = V_USER(122)
  V(24) = V_USER(123)
  V(41) = V_USER(124)
  V(82) = V_USER(125)
  V(70) = V_USER(126)
  V(25) = V_USER(127)
  V(26) = V_USER(128)
  V(27) = V_USER(129)
  V(28) = V_USER(130)
  V(29) = V_USER(131)
  V(30) = V_USER(132)
  V(31) = V_USER(133)
  V(32) = V_USER(134)
  V(33) = V_USER(135)
  V(42) = V_USER(136)
  V(83) = V_USER(137)
  V(34) = V_USER(139)
      
END SUBROUTINE Shuffle_user2kpp

! End of Shuffle_user2kpp function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Shuffle_kpp2user - function to restore concentrations from KPP to USER
!   Arguments :
!      V         - Concentrations of variable species (local)
!      V_USER    - Concentration of variable species in USER's order
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Shuffle_kpp2user ( V, V_USER )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! V_USER - Concentration of variable species in USER's order
  REAL(kind=dp) :: V_USER(NVAR)

  V_USER(1) = V(94)
  V_USER(2) = V(88)
  V_USER(3) = V(123)
  V_USER(4) = V(46)
  V_USER(5) = V(114)
  V_USER(6) = V(106)
  V_USER(7) = V(68)
  V_USER(8) = V(55)
  V_USER(9) = V(116)
  V_USER(10) = V(95)
  V_USER(11) = V(1)
  V_USER(12) = V(2)
  V_USER(13) = V(3)
  V_USER(14) = V(4)
  V_USER(15) = V(5)
  V_USER(16) = V(6)
  V_USER(17) = V(7)
  V_USER(18) = V(8)
  V_USER(19) = V(9)
  V_USER(20) = V(10)
  V_USER(21) = V(11)
  V_USER(22) = V(12)
  V_USER(23) = V(13)
  V_USER(24) = V(14)
  V_USER(25) = V(15)
  V_USER(26) = V(16)
  V_USER(27) = V(111)
  V_USER(28) = V(56)
  V_USER(29) = V(96)
  V_USER(30) = V(102)
  V_USER(31) = V(47)
  V_USER(32) = V(120)
  V_USER(33) = V(49)
  V_USER(34) = V(93)
  V_USER(35) = V(57)
  V_USER(36) = V(135)
  V_USER(37) = V(73)
  V_USER(38) = V(122)
  V_USER(39) = V(76)
  V_USER(40) = V(115)
  V_USER(41) = V(104)
  V_USER(42) = V(92)
  V_USER(43) = V(89)
  V_USER(44) = V(77)
  V_USER(45) = V(79)
  V_USER(46) = V(110)
  V_USER(47) = V(124)
  V_USER(48) = V(105)
  V_USER(49) = V(121)
  V_USER(50) = V(61)
  V_USER(51) = V(99)
  V_USER(52) = V(48)
  V_USER(53) = V(129)
  V_USER(54) = V(131)
  V_USER(55) = V(109)
  V_USER(56) = V(130)
  V_USER(57) = V(54)
  V_USER(58) = V(100)
  V_USER(59) = V(63)
  V_USER(60) = V(117)
  V_USER(61) = V(43)
  V_USER(62) = V(136)
  V_USER(63) = V(139)
  V_USER(64) = V(134)
  V_USER(65) = V(137)
  V_USER(66) = V(132)
  V_USER(67) = V(45)
  V_USER(68) = V(81)
  V_USER(69) = V(108)
  V_USER(70) = V(67)
  V_USER(71) = V(39)
  V_USER(72) = V(107)
  V_USER(73) = V(80)
  V_USER(74) = V(64)
  V_USER(75) = V(113)
  V_USER(76) = V(133)
  V_USER(77) = V(127)
  V_USER(78) = V(69)
  V_USER(79) = V(51)
  V_USER(80) = V(52)
  V_USER(81) = V(128)
  V_USER(82) = V(112)
  V_USER(83) = V(119)
  V_USER(84) = V(72)
  V_USER(85) = V(85)
  V_USER(86) = V(62)
  V_USER(87) = V(103)
  V_USER(88) = V(71)
  V_USER(89) = V(53)
  V_USER(90) = V(44)
  V_USER(91) = V(17)
  V_USER(92) = V(18)
  V_USER(93) = V(19)
  V_USER(94) = V(59)
  V_USER(95) = V(138)
  V_USER(96) = V(90)
  V_USER(97) = V(58)
  V_USER(98) = V(60)
  V_USER(99) = V(40)
  V_USER(100) = V(65)
  V_USER(101) = V(36)
  V_USER(102) = V(37)
  V_USER(103) = V(38)
  V_USER(104) = V(20)
  V_USER(105) = V(21)
  V_USER(106) = V(22)
  V_USER(107) = V(23)
  V_USER(108) = V(50)
  V_USER(109) = V(98)
  V_USER(110) = V(78)
  V_USER(111) = V(101)
  V_USER(112) = V(86)
  V_USER(113) = V(118)
  V_USER(114) = V(74)
  V_USER(115) = V(75)
  V_USER(116) = V(84)
  V_USER(117) = V(125)
  V_USER(118) = V(126)
  V_USER(119) = V(91)
  V_USER(120) = V(66)
  V_USER(121) = V(97)
  V_USER(122) = V(87)
  V_USER(123) = V(24)
  V_USER(124) = V(41)
  V_USER(125) = V(82)
  V_USER(126) = V(70)
  V_USER(127) = V(25)
  V_USER(128) = V(26)
  V_USER(129) = V(27)
  V_USER(130) = V(28)
  V_USER(131) = V(29)
  V_USER(132) = V(30)
  V_USER(133) = V(31)
  V_USER(134) = V(32)
  V_USER(135) = V(33)
  V_USER(136) = V(42)
  V_USER(137) = V(83)
  V_USER(139) = V(34)
      
END SUBROUTINE Shuffle_kpp2user

! End of Shuffle_kpp2user function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! GetMass - compute total mass of selected atoms
!   Arguments :
!      CL        - Concentration of all species (local)
!      Mass      - value of mass balance
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE GetMass ( CL, Mass )

! CL - Concentration of all species (local)
  REAL(kind=dp) :: CL(NSPEC)
! Mass - value of mass balance
  REAL(kind=dp) :: Mass(1)

      
END SUBROUTINE GetMass

! End of GetMass function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE gckpp_Util

