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
      integer(kind=int64)::nCommands,i,j,nAtoms,nAt3,nFrozenAtoms,  &
        nRot,nConstraints,nVib,iCurrentProjector
      integer(kind=int64),dimension(:),allocatable::frozenAtoms
      real(kind=real64),dimension(:),allocatable::cartesians,  &
        atomicMasses,hEVals,tmpVec
      real(kind=real64),dimension(:,:),allocatable::hMat,hEVecs,hMatMW,  &
        hMatProjectorVectors,hMatProjector,RotVecs,tmpMat1,tmpMat2,  &
        tmpMat3
      character(len=512)::matrixFilename,tmpString
      logical::constrainTranslation,constrainRotation
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile
      type(MQC_Variable)::forceConstants
!
!     Format Statements
!
 1000 Format(1x,'Enter VibAnalysis.')
 1010 Format(3x,'Matrix File: ',A,/)
 1100 Format(1x,'nAtoms=',I4)
 1200 Format(1x,'Freezing atom number ',I4)
 2100 Format(1x,'number of rotational DOFs: ',I1)
!
!
      write(IOut,1000)
!
!     Open the Gaussian matrix file and load the number of atomic centers.

      nCommands = command_argument_count()
      if(nCommands.eq.0)  &
        call mqc_error('No command line arguments provided. The input Gaussian matrix file name is required.')
      call get_command_argument(1,matrixFilename)
      call GMatrixFile%load(matrixFilename)
      write(IOut,1010) TRIM(matrixFilename)
      nAtoms = GMatrixFile%getVal('nAtoms')
      write(IOut,1100) nAtoms
!
!     Figure out nAt3, then allocate memory for key arrays.
!
      nAt3 = 3*nAtoms
      Allocate(frozenAtoms(nAtoms),hEVecs(nAt3,nAt3),hEVals(nAt3),  &
        hMatMW(NAt3,NAt3))
!
!     Set frozenAtoms to 0's and 1's. A value of frozenAtoms(i)=1 indicates atom
!     i is frozen.
!
      frozenAtoms  = 0
      nFrozenAtoms = 0
      do i = 2,nCommands
        call get_command_argument(i,tmpString)
        read(tmpString,'(I)') j
        if(j.le.0.or.j.gt.nAtoms) call mqc_error('Invalid atom number given in frozen atom list.')
        frozenAtoms(j) = 1
        nFrozenAtoms = nFrozenAtoms+1
      endDo
      if(extraPrint.or.nFrozenAtoms.gt.0)  &
        call mqc_print(iOut,frozenAtoms,header='frozenAtoms list')
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
!     Determine the number of constraints, allocate space for the projector
!     vectors, and initialize them.
!
      select case(nFrozenAtoms)
      case(0)
        nConstraints = 6
        constrainTranslation = .True.
        constrainRotation = .True.
      case(1)
        nConstraints = nFrozenAtoms*3 + 3
        constrainTranslation = .False.
        constrainRotation = .True.
      case(2:)
        nConstraints = nFrozenAtoms*3
        constrainTranslation = .False.
        constrainRotation = .False.
      case default
        call mqc_error('Determination of constraints is confused.')
      end select
      Allocate(hMatProjectorVectors(nAt3,nConstraints))
      hMatProjectorVectors = float(0)
      iCurrentProjector = 1
!
!     Now, build projectors to remove overall translational degrees of freedom.
!
      if(constrainTranslation) then
        if(Allocated(tmpVec)) deAllocate(tmpVec)
        Allocate(tmpVec(nAt3))
        do i = 1,3
          tmpVec = float(0)
          do j = 0,nAtoms-1
            tmpVec(3*j+i) = float(1)
          endDo
          call massWeighVector(.true.,atomicMasses,tmpVec)
          call mqc_normalizeVector(tmpVec)
          if(iCurrentProjector.gt.nConstraints)  &
            call mqc_error('Logic Error: vibAnalysis attempted to add more constraints than expected.')
          hMatProjectorVectors(:,iCurrentProjector) = tmpVec
          iCurrentProjector = iCurrentProjector+1
        endDo
        deAllocate(tmpVec)
        if(extraPrint)  &
          call mqc_print(IOut,hMatProjectorVectors(:,1:3),header='MW translational projection vector.')
      endIf
!
!     Determine the moments of inertia and principle axes of rotation. Then,
!     build the rotational constaint vectors. As appropriate to the job, add
!     them to hMatProjectorVectors.
!
      if(constrainRotation) then
        Allocate(RotVecs(nAt3,3))
        call momentsOfInertia(iOut,nAtoms,cartesians,atomicMasses,nRot,RotVecs)
        call mqc_normalizeVector(RotVecs(:,1))
        call mqc_normalizeVector(RotVecs(:,2))
        call mqc_normalizeVector(RotVecs(:,3))
        if(extraPrint) then
          write(iOut,2100) nRot
          call mqc_print(iOut,RotVecs,header='(Normalized) Rotational Constraint Vectors')
        endIf
        hMatProjectorVectors(:,iCurrentProjector) = RotVecs(:,1)
        iCurrentProjector = iCurrentProjector+1
        hMatProjectorVectors(:,iCurrentProjector) = RotVecs(:,2)
        iCurrentProjector = iCurrentProjector+1
        hMatProjectorVectors(:,iCurrentProjector) = RotVecs(:,3)
        iCurrentProjector = iCurrentProjector+1
      endIf
!
!     Build the projector based on hMatProjectorVectors.
!
      Allocate(hMatProjector(nAt3,nAt3))
      call mqc_print(iOut,hMatProjectorVectors,header='Projection Vectors')
      hMatProjector = MatMul(hMatProjectorVectors,Transpose(hMatProjectorVectors))
      hMatProjector = unitMatrix(nAt3) - hMatProjector
      call mqc_print(iOut,hMatProjector,header='Projection Matrix -- FINAL')
!
!     Apply the projector to the MW Hessian and diagonalize again.
!
      hMatMW = hMat
      call massWeighMatrix(.false.,atomicMasses,hMatMW)
      hMatMW = MatMul(MatMul(hMatProjector,hMatMW),hMatProjector)
      hMatMW = hMatMW*scaleHess
      if(extraPrint)  &
        call mqc_print(IOut,hMatMW,header='Projected MW Hessian after scaleHess')
      call mySVD(iOut,nAt3,hMatMW,hEVals,hEVecs)
      i = 3+nRot
      if(extraPrint) then
        call mqc_print(iOut,MatMul(TRANSPOSE(hMatProjectorVectors),hEVecs(:,1:i)),header='Overlaps of v and first eVecs')
        call mqc_print(iOut,MatMul(TRANSPOSE(hMatProjectorVectors),hMatProjectorVectors),header='Overlaps of v and v')
        call mqc_print(iOut,MatMul(TRANSPOSE(hMatProjectorVectors),hEVecs(:,i+1:)),header='Overlaps of v and normal modes')
        call mqc_print(IOut,hEVecs,header='EVecs after hMatMW SVD')
        call mqc_print(IOut,hEVals,header='EVals after hMatMW SVD')
      endIf
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
        call mqc_normalizeVector(hEVecs(:,i))
      endDo
      call mqc_print(iOut,hEVals,header='Eigenvalues (cm-1)')
      call mqc_print(IOut,hEVecs,header='Displacements')
!
  999 Continue
      write(iOut,*)' END OF VIBANALYSIS'
      end program VibAnalysis
