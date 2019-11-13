INCLUDE 'vibAnalysisMod.f03'

      Program unitTest
!
!     This program is a unit test program for various sub-programs being written
!     for the vibrational analysis program I'm testing.
!
!     -H. P. Hratchian, 2019.
!
!
!     USE Connections
!
      use VibAnalysisMod
      use mqc_general
      use mqc_gaussian
      use mqc_algebra2
      use iso_fortran_env
!
!     Variable Declarations
!
      implicit none
      integer(kind=int64),parameter::IOut=6
      integer(kind=int64)::nDim
      real(kind=real64),dimension(:),allocatable::eVals
      real(kind=real64),dimension(:,:),allocatable::mat,eVecs
!
!     Format Statements
!
 1000 Format(1x,'Enter z.')
 1010 Format(3x,'Matrix File: ',A,/)
 1100 Format(1x,'nAtoms=',I4)
!
!
      write(IOut,1000)
      nDim = 3
      Allocate(mat(nDim,nDim),eVals(nDim),eVecs(nDim,nDim))
      mat(1,1) =  1
      mat(2,1) = -3
      mat(3,1) =  3
      mat(2,2) =  1
      mat(3,2) =  2
      mat(3,3) = -4
      mat(1,2) = mat(2,1)
      mat(1,3) = mat(3,1)
      mat(2,3) = mat(3,2)
      call mqc_print(iOut,mat,header='mat')
      call mySVD(iOut,nDim,mat,eVals,eVecs)
      call mqc_print(iOut,mat,header='AFTER mySVD, mat')
      call mqc_print(iOut,eVals,header='EigenValues')
      call mqc_print(iOut,eVecs,header='EigenVectors')
!
  999 Continue
      write(iOut,*)' END OF unitTest'
      end program unitTest
