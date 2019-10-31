      Module VibAnalysisMod
      use mqc_general
      use iso_fortran_env
!
      implicit none
      real(kind=real64),parameter::bohr2ang=0.529177d0
      CONTAINS
!
!
!     PROCEDURE OuterProduct
!
      function outerProduct(v1,v2) result(outMat)
      implicit none
      real(kind=real64),dimension(:),intent(in)::v1,v2
      real(kind=real64),dimension(:,:),allocatable,intent(out)::outMat
      integer(kind=int64)::i,j,nDim1,nDim2
!
      nDim1 = Size(v1)
      nDim2 = Size(v2)
!
      Allocate(outMat(nDim1,nDim2))
      do i = 1,nDim1
        do j = 1,nDim2
          outMat(i,j) = v1(i)*v2(j)
        endDo
      endDo
!
      return
      end function outerProduct

!
!     PROCEDURE unitMatrix
!
      function unitMatrix(nDim) result(outMat)
      implicit none
      integer(kind=int64),intent(in)::nDim
      real(kind=real64),dimension(:,:),allocatable,intent(out)::outMat
      integer(kind=int64)::i
!
      Allocate(outMat(nDim,nDim))
      outMat = float(0)
      do i = 1,nDim
        outMat(i,i) = float(1)
      endDo
!
      return
      end function unitMatrix

!
!     PROCEDURE MassWeighVector
!
      subroutine massWeighVector(DoMultiplyOrDivide,AtomicMasses,Vector)
      implicit none
      logical,intent(in)::DoMultiplyOrDivide
      real(kind=real64),dimension(:),intent(in)::AtomicMasses
      real(kind=real64),dimension(:),intent(inOut)::Vector
      integer(kind=int64)::i,j,nAtoms
!
      if(3*Size(AtomicMasses).ne.Size(Vector))  & 
        call mqc_error('MassWeighVector has Vector and AtomicMasses with different lengths.')
      nAtoms = Size(AtomicMasses)
      do i = 1,nAtoms
        j = 3*(i-1)+1
        if(DoMultiplyOrDivide) then
          Vector(j)   = Vector(j)*SQRT(AtomicMasses(i))
          Vector(j+1) = Vector(j+1)*SQRT(AtomicMasses(i))
          Vector(j+2) = Vector(j+2)*SQRT(AtomicMasses(i))
        else
          Vector(j)   = Vector(j)/SQRT(AtomicMasses(i))
          Vector(j+1) = Vector(j+1)/SQRT(AtomicMasses(i))
          Vector(j+2) = Vector(j+2)/SQRT(AtomicMasses(i))
        endIf
      endDo
!
      end subroutine massWeighVector

!
!     PROCEDURE MassWeighMatrix
!
      subroutine massWeighMatrix(DoMultiplyOrDivide,AtomicMasses,Matrix)
      implicit none
      logical,intent(in)::DoMultiplyOrDivide
      real(kind=real64),dimension(:),intent(in)::AtomicMasses
      real(kind=real64),dimension(:,:),intent(inOut)::Matrix
      integer(kind=int64)::i,j,k,nAtoms
!
      if(3*Size(AtomicMasses).ne.Size(Matrix,1))  & 
        call mqc_error('MassWeighMatrix has Matrix and AtomicMasses with different lengths.')
      nAtoms = Size(AtomicMasses)
      do i = 1,nAtoms
        k = 3*(i-1)+1
        if(DoMultiplyOrDivide) then
          Matrix(k,:)   = Matrix(k,:)*SQRT(AtomicMasses(i))
          Matrix(k+1,:) = Matrix(k+1,:)*SQRT(AtomicMasses(i))
          Matrix(k+2,:) = Matrix(k+2,:)*SQRT(AtomicMasses(i))
        else
          Matrix(k,:)   = Matrix(k,:)/SQRT(AtomicMasses(i))
          Matrix(k+1,:) = Matrix(k+1,:)/SQRT(AtomicMasses(i))
          Matrix(k+2,:) = Matrix(k+2,:)/SQRT(AtomicMasses(i))
        endIf
        do j = 1,NAtoms
          k = 3*(j-1)+1
          if(DoMultiplyOrDivide) then
            Matrix(:,k)   = Matrix(:,k)*SQRT(AtomicMasses(i))
            Matrix(:,k+1) = Matrix(:,k+1)*SQRT(AtomicMasses(i))
            Matrix(:,k+2) = Matrix(:,k+2)*SQRT(AtomicMasses(i))
          else
            Matrix(:,k)   = Matrix(:,k)/SQRT(AtomicMasses(i))
            Matrix(:,k+1) = Matrix(:,k+1)/SQRT(AtomicMasses(i))
            Matrix(:,k+2) = Matrix(:,k+2)/SQRT(AtomicMasses(i))
          endIf
        endDo
      endDo
!
      end subroutine massWeighMatrix

!
!     PROCEDURE mySVD
!
      subroutine mySVD(iOut,n,A,eVals,eVecs)
      implicit none
      integer(kind=int64),intent(in)::iOut,n
      real(kind=real64),dimension(n,n),intent(in)::A
      real(kind=real64),dimension(n),intent(out)::eVals
      real(kind=real64),dimension(n,n),intent(out)::eVecs
!
      integer(kind=int64)::iError,lenWork
      real(kind=real64),dimension(:),allocatable::work
      real(kind=real64),dimension(n,n)::copyA,rightEVecs
!
!
      copyA = A
      Allocate(work(1))
      lenWork = -1
      iError = 0
      Call DGESVD('A','A',n,n,copyA,n,eVals,eVecs,n,rightEVecs,n,work,  &
        lenWork,iError)
      lenWork = INT(work(1))
      DeAllocate(work)
      write(IOut,*)' Hrant - lenWork is calculated to be ',lenWork
      Allocate(work(lenWork))
      Call DGESVD('A','A',n,n,copyA,n,eVals,eVecs,n,rightEVecs,n,work,  &
        lenWork,iError)
      DeAllocate(work)
      write(IOut,*)' Hrant - iError = ',iError
      call mqc_print(IOut,rightEVecs,header='mySVD: Right Eigenvectors')
      end subroutine mySVD


!     PROCEDURE momentsOfInertia
!
      subroutine momentsOfInertia(iOut,nAtoms,cartesians,atomicMasses,  &
        nRot,RotVecs)
      implicit none
      integer(kind=int64),intent(in)::iOut,nAtoms
      real(kind=real64),dimension(:),intent(in)::cartesians,atomicMasses
      integer(kind=int64),intent(out)::nRot
      real(kind=real64),dimension(3*nAtoms,3),intent(inOut)::RotVecs
!
      integer(kind=int64)::i,j,iCartOff
      real(kind=real64)::totalMass,cartXP,cartYP,cartZP
      real(kind=real64),dimension(3)::centerOfMass
      real(kind=real64),dimension(:),allocatable::cartesiansCOM
      real(kind=real64),dimension(3)::inertiaEVals
      real(kind=real64),dimension(3,3)::inertiaMat,inertiaEVecs
      real(kind=real64),parameter::Small=1.0d-6
!
!
!     Begin by finding the center of mass.
!
      call mqc_print(iOut,cartesians*bohr2ang,header='input cartesians (A)')
      totalMass = SUM(atomicMasses)
      centerOfMass = float(0)
      do i = 1,nAtoms
        iCartOff = (i-1)*3
        centerOfMass(1) = centerOfMass(1) + atomicMasses(i)*cartesians(iCartOff+1)
        centerOfMass(2) = centerOfMass(2) + atomicMasses(i)*cartesians(iCartOff+2)
        centerOfMass(3) = centerOfMass(3) + atomicMasses(i)*cartesians(iCartOff+3)
      endDo
      centerOfMass = centerOfMass/totalMass
      Allocate(cartesiansCOM(3*nAtoms))
      cartesiansCOM = float(0)
      call mqc_print(iout,cartesians,header='before moving to COM, carts:')
      do i = 1,nAtoms
        iCartOff = (i-1)*3
        write(iout,*)' Hrant - i       =',i
        write(IOut,*)'         iCartOff=',iCartOff
        cartesiansCOM(iCartOff+1) = cartesians(iCartOff+1) - centerOfMass(1)
        cartesiansCOM(iCartOff+2) = cartesians(iCartOff+2) - centerOfMass(2)
        cartesiansCOM(iCartOff+3) = cartesians(iCartOff+3) - centerOfMass(3)
      endDo
      call mqc_print(iOut,centerOfMass*bohr2ang,header='COM (A)')
      call mqc_print(iOut,cartesiansCOM*bohr2ang,header='cartesiansCOM (A)')
!
!     Build the inertia matrix. Then, diagaonalize the matrix to solve for the
!     principal moments and associated eigenvectors.
!
      inertiaMat = float(0)
      do i = 1,nAtoms
        iCartOff = (i-1)*3
        inertiaMat(1,1) = inertiaMat(1,1) + atomicMasses(i)*  &
          (cartesiansCOM(iCartOff+2)**2+cartesiansCOM(iCartOff+3)**2)
        inertiaMat(2,1) = inertiaMat(2,1) - atomicMasses(i)*  &
          cartesiansCOM(iCartOff+1)*cartesiansCOM(iCartOff+2)
        inertiaMat(3,1) = inertiaMat(3,1) - atomicMasses(i)*  &
          cartesiansCOM(iCartOff+1)*cartesiansCOM(iCartOff+3)
        inertiaMat(2,2) = inertiaMat(2,2) + atomicMasses(i)*  &
          (cartesiansCOM(iCartOff+1)**2+cartesiansCOM(iCartOff+3)**2)
        inertiaMat(3,2) = inertiaMat(3,2) - atomicMasses(i)*  &
          cartesiansCOM(iCartOff+2)*cartesiansCOM(iCartOff+3)
        inertiaMat(3,3) = inertiaMat(3,3) + atomicMasses(i)*  &
          (cartesiansCOM(iCartOff+1)**2+cartesiansCOM(iCartOff+2)**2)
      endDo
      inertiaMat(1,2) = inertiaMat(2,1)
      inertiaMat(1,3) = inertiaMat(3,1)
      inertiaMat(2,3) = inertiaMat(3,2)
      call mqc_print(iout,inertiaMat,header='Inertia Matrix:')
      call mySVD(iOut,3,inertiaMat,inertiaEVals,inertiaEVecs)
      call mqc_print(IOut,inertiaEVals,header='inertia Eigenvalues')
      call mqc_print(IOut,inertiaEVecs,header='inertia Eigenvectors')
      do i = 1,3
        if(ABS(inertiaEVals(i)).gt.Small) nRot = nRot+1
      endDo
!
!     Build the three (or two) overall-rotational vectors and determine nRot.
!
      nRot = 0
      do i = 1,3
        if(ABS(inertiaEVals(i)).gt.Small) nRot = nRot+1
      endDo
      if(nRot.ne.2.and.nRot.ne.3)  &
        call mqc_error('Incorrect number of rotational DOF.')
      do j = 1,nAtoms
        iCartOff = (j-1)*3
        cartXP = cartesiansCOM(iCartOff+1)*inertiaEVecs(1,1) +  &
          cartesiansCOM(iCartOff+2)*inertiaEVecs(2,1) +  &
          cartesiansCOM(iCartOff+3)*inertiaEVecs(3,1)
        cartYP = cartesiansCOM(iCartOff+1)*inertiaEVecs(1,2) +  &
          cartesiansCOM(iCartOff+2)*inertiaEVecs(2,2) +  &
          cartesiansCOM(iCartOff+3)*inertiaEVecs(2,3)
        cartZP = cartesiansCOM(iCartOff+1)*inertiaEVecs(1,3) +  &
          cartesiansCOM(iCartOff+2)*inertiaEVecs(2,3) +  &
          cartesiansCOM(iCartOff+3)*inertiaEVecs(3,3)
        RotVecs(iCartOff+1,1) = cartYP*inertiaEVecs(1,3)-cartZP*inertiaEVecs(1,2)
        RotVecs(iCartOff+2,1) = cartYP*inertiaEVecs(2,3)-cartZP*inertiaEVecs(2,2)
        RotVecs(iCartOff+3,1) = cartYP*inertiaEVecs(3,3)-cartZP*inertiaEVecs(3,2)

        RotVecs(iCartOff+1,2) = cartZP*inertiaEVecs(1,1)-cartXP*inertiaEVecs(1,3)
        RotVecs(iCartOff+2,2) = cartZP*inertiaEVecs(2,1)-cartXP*inertiaEVecs(2,3)
        RotVecs(iCartOff+3,2) = cartZP*inertiaEVecs(3,1)-cartXP*inertiaEVecs(3,3)
        if(nRot.eq.3) then
          RotVecs(iCartOff+1,3) = cartXP*inertiaEVecs(1,2)-cartYP*inertiaEVecs(1,1)
          RotVecs(iCartOff+2,3) = cartXP*inertiaEVecs(2,2)-cartYP*inertiaEVecs(2,1)
          RotVecs(iCartOff+3,3) = cartXP*inertiaEVecs(3,2)-cartYP*inertiaEVecs(3,1)
        else
          RotVecs(iCartOff+1,3) = float(0)
          RotVecs(iCartOff+2,3) = float(0)
          RotVecs(iCartOff+3,3) = float(0)
        endIf
      endDo
      call massWeighVector(.true.,atomicMasses,RotVecs(:,1))
      call massWeighVector(.true.,atomicMasses,RotVecs(:,2))
      call massWeighVector(.true.,atomicMasses,RotVecs(:,3))
      call mqc_print(iOut,RotVecs,header='Overall rotation vector:')
!
      end subroutine momentsOfInertia


      End Module VibAnalysisMod




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
      integer(kind=int64)::i,j,nAtoms,nAt3,nRot
      real(kind=real64)::fudge
      real(kind=real64),dimension(:),allocatable::cartesians,  &
        atomicMasses,hEVals,tmpVec
      real(kind=real64),dimension(:,:),allocatable::hMat,hEVecs,hMatMW,  &
        hMatProjector,RotVecs
      character(len=512)::matrixFilename
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile
      type(MQC_Variable)::forceConstants
!
!     Format Statements
!
 1000 Format(1x,'Enter VibAnalysis.')
 1010 Format(3x,'Matrix File: ',A,/)
 1100 Format(1x,'nAtoms=',I4)
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
!     Load the nuclear force constant matrix.
!
      call GMatrixFile%getArray('NUCLEAR FORCE CONSTANTS',mqcVarOut=forceConstants)
      call forceConstants%print(header='force constant matrix')
      hMat = forceConstants
      call mqc_print(IOut,hMat,header='Hessian-FULL')
!
!     Diagonalize the hessian...
!
      call mySVD(iOut,nAt3,hMat,hEVals,hEVecs)
      call mqc_print(IOut,hEVals,header='Eigenvalues')
      call mqc_print(IOut,hEVecs,header='Left Eigenvectors')
!
!     Get the atomic masses then mass-weigh the hessian..
!
      Allocate(atomicMasses(nAtoms))
      atomicMasses = GMatrixFile%getAtomWeights()
      call mqc_print(iout,atomicMasses,header='Atomic Masses')
      hMatMW = hMat
      call massWeighMatrix(.false.,atomicMasses,hMatMW)
!
!     Diagonalize the mass-weighted hessian...
!
      call mySVD(iOut,nAt3,hMatMW,hEVals,hEVecs)
      call mqc_print(IOut,hEVals,header='MW Eigenvalues')
      call mqc_print(IOut,hEVecs,header='MW Left Eigenvectors')
!
!     Un-mass-weigh the eigenvectors and re-print them.
!
      do i = 1,NAt3
        call massWeighVector(.false.,atomicMasses,hEVecs(:,i))
      endDo
      call mqc_print(IOut,hEVecs,header='Un-MW Left Eigenvectors')
!
!     Now, let's build up a projector to remove possible contamination from
!     overall translational degrees of freedom.
!
      Allocate(hMatProjector(nAt3,nAt3))
      if(Allocated(tmpVec)) deAllocate(tmpVec)
      Allocate(tmpVec(nAt3))
      hMatProjector = float(0)
      call mqc_print(IOut,hMatProjector,header='Initial Hessian Projector')
      do i = 1,3
        tmpVec = float(0)
        do j = 0,nAtoms-1
          tmpVec(3*j+i) = float(1)
        endDo
        call massWeighVector(.true.,atomicMasses,tmpVec)
        call mqc_print(IOut,tmpVec,header='MW translational projection vector.')
        hMatProjector = hMatProjector + outerProduct(tmpVec,tmpVec)
        call mqc_print(IOut,hMatProjector,header='Current Hessian Projector')
      endDo
!
!     Determine the moments of inertia and principle axes of rotation.
!
      write(IOut,*)
      write(iOut,*)
      write(iOut,*)' Calling MoI Routine...'
      Allocate(RotVecs(nAt3,3))
      call momentsOfInertia(iOut,nAtoms,cartesians,atomicMasses,nRot,RotVecs)
      write(iOUt,*)' After momentsOfIneria: nRot=',nRot
      call mqc_print(iOut,RotVecs,header='RotVecs in Main Program Unit-1:')
      call mqc_normalizeVector(RotVecs(:,1))
      call mqc_normalizeVector(RotVecs(:,2))
      call mqc_normalizeVector(RotVecs(:,3))
      call mqc_print(iOut,RotVecs,header='RotVecs in Main Program Unit-2:')
      write(IOut,*)
      write(IOut,*)
      write(IOut,*)
      write(iOut,*)' Back from MoI Routine.'
      write(IOut,*)
      do i = 1,nRot
        hMatProjector = hMatProjector + outerProduct(RotVecs(:,i),RotVecs(:,i))
        call mqc_print(IOut,hMatProjector,header='Current Hessian Projector')
      endDo
!
!     Prepare the operator to project out the subspace in hMatProjector.
!
      hMatProjector = unitMatrix(nAt3) - hMatProjector
      call mqc_print(IOut,hMatProjector,header='Final Hessian Projector')
!
!     Now, apply the projector to the MW Hessian and diagonalize again.
!
      hMatMW = MatMul(hMatProjector,hMatMW)
      call mySVD(iOut,nAt3,hMatMW,hEVals,hEVecs)
      call mqc_print(IOut,hEVals,header='MW Eigenvalues')
      call mqc_print(IOut,hEVecs,header='MW Left Eigenvectors')
      fudge = 24162.35735
      hEVals = hEVals*fudge
      call mqc_print(IOut,hEVals,header='MW Eigenvalues (cm-1)')
!
      end program VibAnalysis
