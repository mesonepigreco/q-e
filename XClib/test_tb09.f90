! Copyright (C) 2026 Quantum ESPRESSO group
! Distributed under the GNU General Public License; see License.
!
! Fixed local inputs evaluated independently from the BR/TB09 equations
! (PRL 102, 226401), using a bracketed solution of the BR equation.
! The c=1.2 cases cover signed Laplacians, a gradient extremum, and both
! Libxc spin layouts. They do not depend on a pseudopotential or a fitted gap.
PROGRAM test_tb09
  USE kind_l, ONLY: DP
  USE xc_lib
#if defined(__MPI)
  USE mpi
#endif
  IMPLICIT NONE
  INTEGER, PARAMETER :: n=5
  REAL(DP), PARAMETER :: tol=1.e-9_dp
  REAL(DP) :: rho(n,2), grho(3,n,2), tau(n,2), lapl(n,2)
  REAL(DP) :: ex(n), ec(n), vx(n,2), v2x(n,2), v3x(n,2)
  REAL(DP) :: vc(n,2), v2c(3,n,2), v3c(n,2), ref(n)
  REAL(DP) :: rho_lda(n,2)
  REAL(DP) :: elx(n), elc(n), vlx(n,2), vlc(n,2), unpol_corr(n)
  CHARACTER(LEN=32) :: arg
  INTEGER :: ns, is
#if defined(__MPI)
  INTEGER :: ierr
  CALL MPI_INIT(ierr)
#endif
  CALL get_command_argument(1,arg)
  DO ns=1,2
    CALL xclib_reset_dft()
    CALL xclib_set_dft_from_name('TB09')
    IF (.NOT.xclib_dft_is('TB09')) STOP 1
    IF (xclib_get_ID('MGGA','EXCH')/=208 .OR. &
        xclib_get_ID('MGGA','CORR')/=0 .OR. &
        xclib_get_ID('LDA','EXCH')/=0 .OR. &
        xclib_get_ID('LDA','CORR')/=4) STOP 2
    IF (xclib_get_dft_short()/='TB09') STOP 3
    CALL xclib_init_libxc(ns,.FALSE.)
    CALL set_libxc_ext_param(5,0,1.2_dp)
    rho=.1_dp/REAL(ns,dp)
    tau=.08_dp/REAL(ns,dp)
    grho=0._dp
    lapl=0._dp
    DO is=1,ns
      grho(1,1:2,is)=.05_dp/REAL(ns,dp)
      grho(1,5,is)=.05_dp/REAL(ns,dp)
      lapl(:,is)=[.03_dp,-.03_dp,0._dp,.03_dp,1._dp]/REAL(ns,dp)
    ENDDO
    IF (arg=='missing-laplacian') THEN
      CALL xc_metagcx(n,ns,2*ns-1,rho(:,1:ns),grho(:,:,1:ns),tau(:,1:ns), &
                      ex,ec,vx(:,1:ns),v2x(:,1:ns),v3x(:,1:ns), &
                      vc(:,1:ns),v2c(1:2*ns-1,:,1:ns),v3c(:,1:ns))
      ! The missing-Laplacian CTest must fail, so success here exposes the bug.
      STOP 0
    ENDIF
    CALL xc_metagcx(n,ns,2*ns-1,rho(:,1:ns),grho(:,:,1:ns),tau(:,1:ns), &
                    ex,ec,vx(:,1:ns),v2x(:,1:ns),v3x(:,1:ns), &
                    vc(:,1:ns),v2c(1:2*ns-1,:,1:ns),v3c(:,1:ns),lapl_=lapl(:,1:ns))
    ref=[-0.398671988252366_dp,-0.375393302979790_dp, &
         -0.382634452174017_dp,-0.394468869326934_dp,-0.727706358110888_dp]
    DO is=1,ns
      IF (MAXVAL(ABS(vx(:,is)-ref))>tol) THEN
        WRITE(*,*) 'TB09 mismatch, spin layout:',ns,' channel:',is
        WRITE(*,*) 'Actual: ',vx(:,is)
        WRITE(*,*) 'Reference: ',ref
        STOP 4
      ENDIF
    ENDDO
    IF (MAXVAL(ABS(v2x(:,1:ns)))>tol .OR. MAXVAL(ABS(v3x(:,1:ns)))>tol) STOP 5
    IF (MAXVAL(ABS(ex))>tol .OR. MAXVAL(ABS(ec))>tol) STOP 6
    IF (MAXVAL(ABS(vc(:,1:ns)))>tol .OR. MAXVAL(ABS(v3c(:,1:ns)))>tol) STOP 7
    ! The LDA correlation is evaluated separately, exactly as in PW.
    ! The LDA interface takes total density/magnetization, not up/down.
    rho_lda(:,1)=SUM(rho(:,1:ns),DIM=2)
    IF (ns==2) rho_lda(:,2)=rho(:,1)-rho(:,2)
    CALL xc(n,ns,ns,rho_lda(:,1:ns),elx,elc,vlx(:,1:ns),vlc(:,1:ns))
    IF (MAXVAL(ABS(elx))>tol .OR. MAXVAL(ABS(vlx(:,1:ns)))>tol) STOP 8
    IF (ns==1) unpol_corr=vlc(:,1)
    DO is=1,ns
      IF (MAXVAL(ABS(vlc(:,is)-unpol_corr))>tol) STOP 9
    ENDDO
    IF (ns==2) THEN
      ! Unequal spin channels expose an incorrect polarized-array layout.
      rho(:,1)=.06_dp; rho(:,2)=.04_dp
      tau(:,1)=.05_dp; tau(:,2)=.03_dp
      grho=0._dp
      grho(1,:,1)=SQRT(.001_dp); grho(1,:,2)=SQRT(.0005_dp)
      lapl(:,1)=.015_dp; lapl(:,2)=-.01_dp
      CALL xc_metagcx(n,2,3,rho,grho,tau,ex,ec,vx,v2x,v3x,vc,v2c,v3c,lapl_=lapl)
      IF (MAXVAL(ABS(vx(:,1)+.446589368326217_dp))>tol .OR. &
          MAXVAL(ABS(vx(:,2)+.326535662094257_dp))>tol) STOP 10
    ENDIF
    CALL xclib_finalize_libxc()
  ENDDO
  WRITE(*,*) 'TB09 local-potential and alias checks passed'
#if defined(__MPI)
  CALL MPI_FINALIZE(ierr)
#endif
END PROGRAM test_tb09
