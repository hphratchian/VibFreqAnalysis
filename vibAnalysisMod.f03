      Module VibAnalysisMod
      use mqc_general
      use iso_fortran_env
!
      implicit none
      integer(kind=int64),parameter::IOut=6
      real(kind=real64),parameter::bohr2ang=0.529177d0,  &
        scale2wavenumber=48169.11381488435d0,  &
        scaleHess=548.5799089940967d0
      logical::extraPrint=.False.
      CONTAINS
!
!
!     PROCEDURE addConstraintVector
!
      subroutine addConstraintVector(iCurrentProjector,nConstraints,  &
        newVector,hMatProjectorVectors)
!
!     This subroutine adds the input vector newVector to the list of projectors
!     stored in hMatProjectorVectors. The dummy argument iCurrentProjector
!     should be initialized to zero before the first call to this routine. The
!     value in iCurrentProjector is then updated here and should not be modified
!     outside this routine. Dummy argument nConstraints is the number of
!     expected total constraints.
!
      implicit none
      integer(kind=int64),intent(inOut)::iCurrentProjector
      integer(kind=int64),intent(in)::nConstraints
      real(kind=real64),dimension(:),intent(in)::newVector
      real(kind=real64),dimension(:,:),intent(inOut)::hMatProjectorVectors
!
!     Do the work...
!
      if(iCurrentProjector.gt.nConstraints)  &
        call mqc_error('Logic Error: vibAnalysis attempted to add more constraints than expected.')
      hMatProjectorVectors(:,iCurrentProjector) = newVector
      iCurrentProjector = iCurrentProjector+1
!
      return
      end subroutine addConstraintVector


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
          call mqc_print(vector,iOut,header='Current Vector')
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
      call mqc_print(vector,iOut,header='FINAL Vector')
      write(6,*)
      write(6,*)
!
      return
      end subroutine orthogonalizeVector

!
!     PROCEDURE vecMatVec
!
      function VecMatVec(vector,matrix) result(outVal)
      implicit none
      real(kind=real64),dimension(:),intent(in)::vector
      real(kind=real64),dimension(:,:),intent(in)::matrix
      real(kind=real64),intent(out)::outVal
      real(kind=real64),dimension(:),allocatable::tmpVector
      integer(kind=int64)::i,j,nDim,nDim1,nDim2
!
      nDim  = Size(vector)
      nDim1 = Size(matrix,1)
      nDim2 = Size(matrix,2)
      if(nDim.ne.nDim1.or.nDim1.ne.nDim2) call mqc_error('VecMatVec: Problem with nDims.')
      Allocate(tmpVector(nDim))
      tmpVector = MatMul(matrix,vector)
      outVal = dot_product(vector,tmpVector)
!
      return
      end function VecMatVec

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
      return
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
      return
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
      integer(kind=int64)::i,iError,lenWork
      real(kind=real64),dimension(:),allocatable::work
      real(kind=real64),dimension(n,n)::copyA,leftEVecs,rightEVecs
!
!
      copyA = A
      Allocate(work(1))
      lenWork = -1
      iError = 0
      Call DGESVD('A','A',n,n,copyA,n,eVals,leftEVecs,n,rightEVecs,n,work,  &
        lenWork,iError)
      lenWork = INT(work(1))
      DeAllocate(work)
      Allocate(work(lenWork))
      Call DGESVD('A','A',n,n,copyA,n,eVals,leftEVecs,n,rightEVecs,n,work,  &
        lenWork,iError)
      DeAllocate(work)
!
!     Resort the eigenvectors...
!
      do i = 0,n-1
        EVecs(:,i+1) = leftEVecs(:,n-i)
      endDo
!
!     Go through all the eigenvectors and determine eigenvalues and fix the
!     phase of each vector so that the largest magnitude element is positive.
!
      Allocate(work(n))
      do i = 1,n
        call mqc_vectorPhase(EVecs(:,i),.TRUE.)
        work = MatMul(A,EVecs(:,i))
        eVals(i) = dot_product(EVecs(:,i),work)
      endDo
      if(extraPrint) then
        call mqc_print(EVecs,iOut,header='mySVD: Left Eigenvectors')
        call mqc_print(rightEVecs,iOut,header='mySVD: Right Eigenvectors')
      endIf
      DeAllocate(work)
!
      return
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
      call mqc_print(cartesians,iOut,header='input cartesians (au)')
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
      if(extraPrint)  &
        call mqc_print(cartesians,iOut,header='before moving to COM, carts:')
      do i = 1,nAtoms
        iCartOff = (i-1)*3
        cartesiansCOM(iCartOff+1) = cartesians(iCartOff+1) - centerOfMass(1)
        cartesiansCOM(iCartOff+2) = cartesians(iCartOff+2) - centerOfMass(2)
        cartesiansCOM(iCartOff+3) = cartesians(iCartOff+3) - centerOfMass(3)
      endDo
      if(extraPrint) then
        call mqc_print(centerOfMass,iOut,header='COM (au)')
        call mqc_print(cartesiansCOM,iOut,header='cartesiansCOM (au)')
      endIf
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
      call mySVD(iOut,3,inertiaMat,inertiaEVals,inertiaEVecs)
      if(extraPrint) then
        call mqc_print(inertiaMat,iOut,header='Inertia Matrix:')
        call mqc_print(inertiaEVals,iOut,header='inertia Eigenvalues')
        call mqc_print(inertiaEVecs,iOut,header='inertia Eigenvectors')
      endIf
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
          cartesiansCOM(iCartOff+3)*inertiaEVecs(3,2)
        cartZP = cartesiansCOM(iCartOff+1)*inertiaEVecs(1,3) +  &
          cartesiansCOM(iCartOff+2)*inertiaEVecs(2,3) +  &
          cartesiansCOM(iCartOff+3)*inertiaEVecs(3,3)
        RotVecs(iCartOff+1,1) = cartYP*inertiaEVecs(1,3)-cartZP*inertiaEVecs(1,2)
        RotVecs(iCartOff+2,1) = cartYP*inertiaEVecs(2,3)-cartZP*inertiaEVecs(2,2)
        RotVecs(iCartOff+3,1) = cartYP*inertiaEVecs(3,3)-cartZP*inertiaEVecs(3,2)
        RotVecs(iCartOff+1,2) = cartZP*inertiaEVecs(1,1)-cartXP*inertiaEVecs(1,3)
        RotVecs(iCartOff+2,2) = cartZP*inertiaEVecs(2,1)-cartXP*inertiaEVecs(2,3)
        RotVecs(iCartOff+3,2) = cartZP*inertiaEVecs(3,1)-cartXP*inertiaEVecs(3,3)
        RotVecs(iCartOff+1,3) = cartXP*inertiaEVecs(1,2)-cartYP*inertiaEVecs(1,1)
        RotVecs(iCartOff+2,3) = cartXP*inertiaEVecs(2,2)-cartYP*inertiaEVecs(2,1)
        RotVecs(iCartOff+3,3) = cartXP*inertiaEVecs(3,2)-cartYP*inertiaEVecs(3,1)
      endDo
      call massWeighVector(.true.,atomicMasses,RotVecs(:,1))
      call massWeighVector(.true.,atomicMasses,RotVecs(:,2))
      call massWeighVector(.true.,atomicMasses,RotVecs(:,3))
      if(extraPrint) call mqc_print(RotVecs,iOut,header='Overall rotation vectors:')
!
      return
      end subroutine momentsOfInertia


      End Module VibAnalysisMod
