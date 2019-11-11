      Module VibAnalysisMod
      use mqc_general
      use iso_fortran_env
!
      implicit none
      real(kind=real64),parameter::bohr2ang=0.529177d0
      CONTAINS
!
!
!     PROCEDURE OrthogonalizeVector
!
      subroutine orthogonalizeVector(nHave,vectorList,vector)
!
!     This subroutine takes a vector and orthogonalizes is relative to the nHave
!     column vectors in vectorList.
!
      implicit none
      integer(kind=int64),intent(in)::nHave
      real(kind=real64),dimension(:,:),intent(in)::vectorList
      real(kind=real64),dimension(:),intent(inOut)::vector
      integer(kind=int64)::i,iTry,nDim
      real(kind=real64)::tmpMagnitude
      real(kind=real64),dimension(:),allocatable::tmpVector
      real(kind=real64),parameter::Small=1.0d-6
!
!     Start by determining nDim. Then do a couple santity checks.
!
      if(nHave.lt.1) call mqc_error('OrthogonalizeVector: nHave < 1.')
      nDim = Size(vector)
      if(nDim.ne.Size(vectorList,1)) call mqc_error('OrthongalizeVector: vectorList has wrong dimension.')
      Allocate(tmpVector(nDim))
!
!     Begin by initializing the input vector with a 1 in a trial location.
      write(6,*)
      write(6,*)
      write(6,*)' OrthogonalizeVector: nHave=',nHave
      do iTry = 1,nDim
        do i = 1,nHave
          call mqc_normalizeVector(vectorList(:,i))
          call mqc_print(6,vector,header='Current Vector')
          tmpMagnitude = dot_product(vectorList(:,i),vector)
          vector = vector - tmpMagnitude*vectorList(:,i)
          tmpMagnitude = dot_product(vector,vector)
          if(tmpMagnitude.lt.Small) exit
        endDo
        if(tmpMagnitude.gt.Small) exit
        vector = tmpVector
        vector(iTry) = vector(iTry) + float(1)
      endDo
      if(tmpMagnitude.le.Small) call mqc_error('Failed making orthog vector.')
      call mqc_normalizeVector(vector)
      call mqc_print(6,vector,header='FINAL Vector')
      write(6,*)
      write(6,*)
!
      return
      end subroutine orthogonalizeVector

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
      integer(kind=int64)::i,j,iAtom,jAtom,nAtoms,nAt3
!
      if(3*Size(AtomicMasses).ne.Size(Matrix,1))  & 
        call mqc_error('MassWeighMatrix has Matrix and AtomicMasses with different lengths.')
      nAtoms = Size(AtomicMasses)
      nAt3 = 3*nAtoms
      do i = 1,nAt3
        iAtom = (i+2)/3
        do j = 1,nAt3
          jAtom = (j+2)/3
          if(DoMultiplyOrDivide) then
            Matrix(i,j) = Matrix(i,j)*SQRT(AtomicMasses(iAtom)*AtomicMasses(jAtom))
          else
            Matrix(i,j) = Matrix(i,j)/SQRT(AtomicMasses(iAtom)*AtomicMasses(jAtom))
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
      write(IOut,*)
      write(IOut,*)
      write(IOut,*)
      write(IOut,*)
      call mqc_print(iOut,cartesians,header='input cartesians (au)')
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
      call mqc_print(iOut,centerOfMass,header='COM (au)')
      call mqc_print(iOut,cartesiansCOM,header='cartesiansCOM (au)')
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
      call mqc_print(iOut,RotVecs,header='Overall rotation vectors:')
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
      call mqc_print(IOut,hMatMW,header='Hessian-FULL after MW''ing')
      hMatMW = hMatMW*scaleHess
      call mqc_print(IOut,hMatMW,header='Hessian-FULL after scaleHess')
!
!     Diagonalize the mass-weighted hessian. This is done prior to projection of
!     translation/rotation/frozen-atom constrains.
!
      call mySVD(iOut,nAt3,hMatMW,hEVals,hEVecs)
      call mqc_print(IOut,hEVals,header='MW Eigenvalues')
      call mqc_print(IOut,hEVecs,header='MW Left Eigenvectors')
      hEVals = hEVals*scale2wavenumber
      hEVals = SIGN(SQRT(ABS(hEVals)),hEVals)
      call mqc_print(IOut,hEVals,header='MW Eigenvalues (cm-1)')
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
      call mqc_print(IOut,hMatProjectorVectors(:,1:3),header='MW translational projection vector.')
!
!     Determine the moments of inertia and principle axes of rotation.
!
      write(IOut,*)
      write(iOut,*)
      write(iOut,*)' Calling MoI Routine...'
      Allocate(RotVecs(nAt3,3))
      call momentsOfInertia(iOut,nAtoms,cartesians,atomicMasses,nRot,RotVecs)
      write(iOUt,*)' After momentsOfIneria: nRot=',nRot
      call mqc_print(iOut,RotVecs,header='RotVecs')
      call mqc_normalizeVector(RotVecs(:,1))
      call mqc_normalizeVector(RotVecs(:,2))
      call mqc_normalizeVector(RotVecs(:,3))
      call mqc_print(iOut,RotVecs,header='RotVecs after normalization')
      hMatProjectorVectors(:,4) = RotVecs(:,1)
      hMatProjectorVectors(:,5) = RotVecs(:,2)
      hMatProjectorVectors(:,6) = RotVecs(:,3)
!
!     Prepare the operator to project out the subspace in hMatProjector.
!
      Allocate(hMatProjector(nAt3,nAt3))
      call mqc_print(iOut,hMatProjectorVectors,header='Projection Vectors')
      hMatProjector = MatMul(hMatProjectorVectors,Transpose(hMatProjectorVectors))
      call mqc_print(iOut,hMatProjector,header='Projection Matrix -- Version 0')
      hMatProjector = unitMatrix(nAt3) - hMatProjector
      call mqc_print(iOut,hMatProjector,header='Projection Matrix -- FINAL')

!hph+
!      write(iOut,*)
!      write(iOut,*)' Hrant - Here''s the first way of building the projector...'
!      hMatProjector = float(0)
!      do i = 1,6
!        hMatProjector = hMatProjector +  &
!          outerProduct(hMatProjectorVectors(:,i),hMatProjectorVectors(:,i))
!!        hMatProjector =   &
!!          outerProduct(hMatProjectorVectors(:,i),hMatProjectorVectors(:,i))
!        write(iOut,*)' Projector ',i,'...'
!        call mqc_print(iOut,hMatProjector,header='     Current Projection Matrix')
!      endDo
!      write(iOut,*)
!      write(iOut,*)
!      write(iOut,*)
!!hph      hMatProjector = unitMatrix(nAt3) - hMatProjector
!      call mqc_print(iOut,hMatProjector,header='FINAL Projection Matrix -- Version 1')
!
!      hMatProjector = MatMul(hMatProjectorVectors,Transpose(hMatProjectorVectors))
!      call mqc_print(iOut,hMatProjector,header='FINAL Projection Matrix -- Version 2')
!hph-

!hph+
!      hMatProjector = unitMatrix(nAt3) -  &
!        MatMul(hMatProjectorVectors(:,1:5),Transpose(hMatProjectorVectors(:,1:5)))
!      call mqc_print(IOut,hMatProjector,header='Final Hessian Projector, hMatProjector')
!      hMatProjector = float(0)
!      do i = 1,6
!        hMatProjector = hMatProjector +  &
!          outerProduct(hMatProjectorVectors(:,1),hMatProjectorVectors(:,1))
!      endDo
!      hMatProjector = unitMatrix(nAt3) - hMatProjector
!      call mqc_print(iOut,hMatProjector,header='hMatProjector, again...')
!hph-

!hph      goto 999
      goto 800

!
!     Build the projection matrix into tmpMat1.
!
      Allocate(tmpMat1(nAt3,nAt3))
      tmpMat1 = float(0)
      tmpMat1(:,1:3+nRot) = hMatProjectorVectors
      write(iOut,*)
      write(iOut,*)
      call mqc_print(iOut,tmpMat1,header='tmpMat1 before Gram-Schmidt...')
      do i = nRot+4,NAt3
        write(IOut,*)
        write(iOUt,*)
        write(IOut,*)' Hrant'
        write(iOut,*)' Orthogonalizing vector i=',i
        tmpMat1(i,i) = float(1)
        if(i.eq.(nRot+4)) then
          tmpMat1(:,i) = float(0)
          tmpMat1(2,i) =  0.256965
          tmpMat1(5,i) = -0.511850
          tmpMat1(6,i) =  0.452762
          tmpMat1(8,i) = -0.511850
          tmpMat1(9,i) = -0.452762
        endIf
        call orthogonalizeVector(i-1,tmpMat1(:,1:i-1),tmpMat1(:,i))
        call mqc_print(iOut,tmpMat1,header='tmpMat1')
      endDo
      write(iOut,*)
      write(iOut,*)
      write(iOut,*)' Check of normalization of vector space...'
      do i = 1,nAt3
        write(iOut,*)' Hrant - i=',i,'   < i | i >=',dot_product(tmpMat1(:,i),tmpMat1(:,i))
      endDo
      call mqc_print(iOut,MatMul(Transpose(tmpMat1),tmpMat1),header='Orthonormality Check')


  800 Continue

!hph+
!!
!!     Project the MW Hessian into the sub-space orthogonal to the constrained
!!     sub-space. Then, diagonalize the new matrix and report eigenvalues.
!!
!      nVib = nAt3-3-nRot
!      Allocate(tmpMat2(nVib,nVib),tmpMat3(nVib,nVib))
!
!      hMatMW = hMat
!      call massWeighMatrix(.false.,atomicMasses,hMatMW)
!      hMatMW = hMatMW*scaleHess
!      call mqc_print(iOut,hMatMW,header='hMatMW before projection...')
!
!      tmpMat2 = MatMul(  &
!        MatMul(Transpose(tmpMat1(:,3+nRot+1:nAt3)),hMatMW),  &
!        tmpMat1(:,3+nRot+1:nAt3))
!
!      tmpMat2 = MatMul(MatMul(hMatProjector,hMatMW),hMatProjector)
!      call mqc_print(iOut,tmpMat2,header='tmpMat2')
!      write(iOut,*)
!
!      hEVals = float(0)
!      call mySVD(iOut,nVib,tmpMat2,hEVals(1:nVib),tmpMat3)
!      call mqc_print(IOut,hEVals,header='EVals after sub-hMatMW SVD')
!      hEVals = hEVals*scale2wavenumber
!      hEVals = SIGN(SQRT(ABS(hEVals)),hEVals)
!      call mqc_print(IOut,hEVals,header='MW Eigenvalues (cm-1)')
!hph-


!
!     Apply the projector to the MW Hessian and diagonalize again.
!
      write(iOut,*)' Hrant - NEW CODE...'
      hMatMW = hMat
      call massWeighMatrix(.false.,atomicMasses,hMatMW)
      hMatMW = MatMul(MatMul(hMatProjector,hMatMW),hMatProjector)
      hMatMW = hMatMW*scaleHess
      call mqc_print(IOut,hMatMW,header='Projected MW Hessian after scaleHess')
      call mySVD(iOut,nAt3,hMatMW,hEVals,hEVecs)
      call mqc_print(IOut,hEVals,header='EVals after hMatMW SVD')
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
