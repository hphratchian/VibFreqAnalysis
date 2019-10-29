      Module VibAnalysisMod
      use mqc_general
      use iso_fortran_env
!
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
      integer(kind=int64)::i,nAtoms,nAt3
      real(kind=real64),dimension(:),allocatable::atomicMasses,hEVals
      real(kind=real64),dimension(:,:),allocatable::hMat,hEVecs,hMatMW
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
!     Now, we need to build projectors to remove possible contamination from
!     overall rotational degrees of freedom.
!




      end program VibAnalysis
