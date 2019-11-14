INCLUDE 'vibAnalysisMod.f03'
      Program VibAnalysis
!
!     This program carries out a vibrational analysis using force constants
!     provided from a Gaussian job (via a matrix file) and atomic masses given
!     by the user at a set of prompts.
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
      integer(kind=int64)::i,j,nAtoms,nAt3,nRot,nVib
      real(kind=real64),parameter::scaleHess=548.5799089940967d0,  &
        scale2wavenumber = 48169.11381488435d0
      real(kind=real64),dimension(:),allocatable::cartesians,  &
        atomicMasses,hEVals,tmpVec
      real(kind=real64),dimension(:,:),allocatable::hMat,hEVecs,hMatMW,  &
        hMatProjectorVectors,hMatProjector,RotVecs,tmpMat1,tmpMat2,  &
        tmpMat3
      character(len=512)::matrixFilename
      logical::extraPrint=.false.
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile
      type(MQC_Variable)::forceConstants
!
!     Format Statements
!
 1000 Format(1x,'Enter VibAnalysis.')
 1010 Format(3x,'Matrix File: ',A,/)
 1100 Format(1x,'nAtoms=',I4)
 2100 Format(1x,'number of rotational DOFs: ',I1)
!
!
      write(IOut,1000)
!
!     Open the Gaussian matrix file and load the number of alpha electrons
!     (nAlpha), number of beta electrons (nBeta), and number of MOs (nBasis).
!
      call get_command_argument(1,matrixFilename)
      call GMatrixFile%load(matrixFilename)
      write(IOut,1010) TRIM(matrixFilename)
      nAtoms = GMatrixFile%getVal('nAtoms')
      write(IOut,1100) nAtoms
!
!     Figure out nAt3, then allocate memory for intrinsic arrays.
!
      nAt3 = 3*nAtoms
      Allocate(hEVecs(nAt3,nAt3),hEVals(nAt3),hMatMW(NAt3,NAt3))
!
!     Load the Cartesian coordinates.
!
      cartesians = GMatrixFile%getAtomCarts()
!
!     Load the nuclear force constant matrix. Note that for numerical stability
!     reasons, we precondition the force constant matrix by the first part of
!     the conversion from AU to wavenumbers.
!
      call GMatrixFile%getArray('NUCLEAR FORCE CONSTANTS',mqcVarOut=forceConstants)
      call forceConstants%print(header='force constant matrix')
      hMat = forceConstants
      call mqc_print(IOut,hMat,header='Hessian-FULL')
!
!     Get the atomic masses then mass-weigh the hessian..
!
      Allocate(atomicMasses(nAtoms))
      atomicMasses = GMatrixFile%getAtomWeights()
      call mqc_print(iout,atomicMasses,header='Atomic Masses')
      hMatMW = hMat
      call massWeighMatrix(.false.,atomicMasses,hMatMW)
      if(extraPrint) call mqc_print(IOut,hMatMW,header='Hessian-FULL after MW''ing')
      hMatMW = hMatMW*scaleHess
      if(extraPrint) call mqc_print(IOut,hMatMW,header='Hessian-FULL after scaleHess')
!
!     Diagonalize the mass-weighted hessian. This is done prior to projection of
!     translation/rotation/frozen-atom constrains.
!
      call mySVD(iOut,nAt3,hMatMW,hEVals,hEVecs)
      hEVals = hEVals*scale2wavenumber
      hEVals = SIGN(SQRT(ABS(hEVals)),hEVals)
      call mqc_print(IOut,hEVals,header='Initial MW Eigenvalues (cm-1)')
!
!     Allocate space for the projector vectors and initialize them.
!
      Allocate(hMatProjectorVectors(nAt3,6))
      hMatProjectorVectors = float(0)
!
!     Now, build up a projector to remove possible contamination from overall
!     translational degrees of freedom.
!
      if(Allocated(tmpVec)) deAllocate(tmpVec)
      Allocate(tmpVec(nAt3))
      do i = 1,3
        tmpVec = float(0)
        do j = 0,nAtoms-1
          tmpVec(3*j+i) = float(1)
        endDo
        call massWeighVector(.true.,atomicMasses,tmpVec)
        call mqc_normalizeVector(tmpVec)
        hMatProjectorVectors(:,i) = tmpVec
      endDo
      if(extraPrint)  &
        call mqc_print(IOut,hMatProjectorVectors(:,1:3),header='MW translational projection vector.')
!
!     Determine the moments of inertia and principle axes of rotation. Then,
!     build the rotational constaint vectors. As appropriate to the job, add
!     them to hMatProjectorVectors.
!
      Allocate(RotVecs(nAt3,3))
      call momentsOfInertia(iOut,nAtoms,cartesians,atomicMasses,nRot,RotVecs)
      call mqc_normalizeVector(RotVecs(:,1))
      call mqc_normalizeVector(RotVecs(:,2))
      call mqc_normalizeVector(RotVecs(:,3))
      if(extraPrint) then
        write(iOut,2100) nRot
        call mqc_print(iOut,RotVecs,header='(Normalized) Rotational Constraint Vectors')
      endIf
      hMatProjectorVectors(:,4) = RotVecs(:,1)
      hMatProjectorVectors(:,5) = RotVecs(:,2)
      hMatProjectorVectors(:,6) = RotVecs(:,3)
!
!     Build the projector based on hMatProjectorVectors.
!
      Allocate(hMatProjector(nAt3,nAt3))
      call mqc_print(iOut,hMatProjectorVectors,header='Projection Vectors')
      hMatProjector = MatMul(hMatProjectorVectors,Transpose(hMatProjectorVectors))
      call mqc_print(iOut,hMatProjector,header='Projection Matrix -- Version 0')
      hMatProjector = unitMatrix(nAt3) - hMatProjector
      call mqc_print(iOut,hMatProjector,header='Projection Matrix -- FINAL')
!
!     Apply the projector to the MW Hessian and diagonalize again.
!
      write(iOut,*)' Hrant - NEW CODE...'
      hMatMW = hMat
      call massWeighMatrix(.false.,atomicMasses,hMatMW)
      hMatMW = MatMul(MatMul(hMatProjector,hMatMW),hMatProjector)
      hMatMW = hMatMW*scaleHess
      call mqc_print(IOut,hMatMW,header='Projected MW Hessian after scaleHess')
      write(iOut,*)
      write(iOut,*)
      write(iOut,*)' Here are my <v|H|v> tests...'
      do i = 1,6
        write(iOut,*)' i = ',i
        write(iOut,*)' <v|H|v> = ',VecMatVec(hMatProjectorVectors(:,i),hMatMW)
      endDo
      write(iOut,*)
      write(iOut,*)
      call mySVD(iOut,nAt3,hMatMW,hEVals,hEVecs)
      write(iOut,*)
      write(iOut,*)
      i = 3+nRot
      call mqc_print(iOut,MatMul(TRANSPOSE(hMatProjectorVectors),hEVecs(:,1:i)),header='Overlaps of v and first eVecs')
      call mqc_print(iOut,MatMul(TRANSPOSE(hMatProjectorVectors),hMatProjectorVectors),header='Overlaps of v and v')
      call mqc_print(iOut,MatMul(TRANSPOSE(hMatProjectorVectors),hEVecs(:,i+1:)),header='Overlaps of v and normal modes')
      write(iOut,*)
      write(iOut,*)
      
      call mqc_print(IOut,hEVecs,header='EVecs after hMatMW SVD')
      call mqc_print(IOut,hEVals,header='EVals after hMatMW SVD')

!hph+
      call mqc_error('STOP 3')
!hph-

      hEVals = hEVals*scale2wavenumber
      hEVals = SIGN(SQRT(ABS(hEVals)),hEVals)
      if(extraPrint) then
        call mqc_print(IOut,hEVals,header='MW Eigenvalues (cm-1)')
        call mqc_print(IOut,hEVecs,header='MW Left Eigenvectors')
      endIf
!
!     Un-mass-weigh the eigenvectors and re-print them. Then, write out the
!     eigenvalues converted to wavenumbers.
!
      do i = 1,NAt3
        call massWeighVector(.false.,atomicMasses,hEVecs(:,i))
      endDo
      call mqc_print(IOut,hEVecs,header='Un-MW Left Eigenvectors')
      call mqc_print(iOut,hEVals,header='Eigenvalues (cm-1)')
!
  999 Continue
      write(iOut,*)' END OF VIBANALYSIS'
      end program VibAnalysis
