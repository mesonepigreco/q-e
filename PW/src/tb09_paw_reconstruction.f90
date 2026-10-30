! Included in paw_onecenter. Dedicated multiplicative-potential PAW transform.
! No ABINIT implementation is used. Frozen core and partial-wave conventions
! are the QE UPF conventions; tau = 1/2 sum |grad(phi)|**2 in Hartree units.
!
! A smooth partition of unity integrates the ORIGINAL nonlinear |grad n|/n
! after reconstructing the density, including the smooth crystal environment.
! The smooth reference potential is arbitrary inside augmentation spheres.
! Its actual FFT interpolation, not a second model evaluated on pseudo atomic
! densities, is subtracted in the one-center matrix elements.

  SUBROUTINE PAW_tb09_reset()
    IF(ALLOCATED(tb09_atoms)) DEALLOCATE(tb09_atoms)
    IF(ALLOCATED(tb09_pw_grid)) DEALLOCATE(tb09_pw_grid)
    IF(ALLOCATED(tb09_neighbors)) DEALLOCATE(tb09_neighbors)
    tb09_coreless = .FALSE.
  END SUBROUTINE

  SUBROUTINE tb09_derivative(r,f,df)
    REAL(DP), INTENT(IN) :: r(:),f(:)
    REAL(DP), INTENT(OUT) :: df(:)
    INTEGER :: i,j,k,lo,hi,n
    REAL(DP) :: w
    n=SIZE(r)
    DO i=1,n
       lo=MAX(1,MIN(i-2,n-4)); hi=lo+4
       df(i)=0._dp
       DO j=lo,hi
          IF(j==i) CYCLE
          w=1._dp/(r(j)-r(i))
          DO k=lo,hi
             IF(k==j.OR.k==i) CYCLE
             w=w*(r(i)-r(k))/(r(j)-r(k))
          ENDDO
          df(i)=df(i)+w*(f(j)-f(i))
       ENDDO
    ENDDO
  END SUBROUTINE

  SUBROUTINE tb09_tau_lm(r,ur,bec,indv,nhtolm,lll,ap,lpx,lpl,tk)
    ! Contraction of positive kinetic density in real spherical harmonics.
    ! bec contains twice the off-diagonal occupations, as in PAW_rho_lm.
    REAL(DP), INTENT(IN) :: r(:),ur(:,:),bec(:),ap(:,:,:)
    INTEGER, INTENT(IN) :: indv(:),nhtolm(:),lll(:),lpx(:,:),lpl(:,:,:)
    REAL(DP), INTENT(OUT) :: tk(:,:)
    REAL(DP) :: dur(SIZE(ur,1),SIZE(ur,2)),angular,pref
    INTEGER :: nb,mb,ih,jh,ij,li,lj,lp,lm,ll
    DO nb=1,SIZE(ur,2)
       CALL tb09_derivative(r,ur(:,nb),dur(:,nb))
    ENDDO
    tk=0._dp; ij=0
    DO ih=1,SIZE(indv)
       nb=indv(ih); li=lll(nb)
       DO jh=ih,SIZE(indv)
          ij=ij+1; mb=indv(jh); lj=lll(mb)
          DO lp=1,lpx(nhtolm(jh),nhtolm(ih))
             lm=lpl(nhtolm(jh),nhtolm(ih),lp); ll=INT(SQRT(REAL(lm-1,dp)))
             angular=0.5_dp*REAL(li*(li+1)+lj*(lj+1)-ll*(ll+1),dp)
             pref=0.5_dp*bec(ij)*ap(lm,nhtolm(ih),nhtolm(jh))
             tk(:,lm)=tk(:,lm)+pref*(dur(:,nb)*dur(:,mb)+angular*ur(:,nb)*ur(:,mb)/r**2)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE

  SUBROUTINE tb09_radial(r,f,df,x,v,dv,d2v)
    REAL(DP), INTENT(IN) :: r(:),f(:),df(:),x
    REAL(DP), INTENT(OUT) :: v,dv,d2v
    INTEGER :: lo,hi,mid
    REAL(DP) :: h,t,a,b,c,d
    lo=1; hi=SIZE(r)
    IF(x>r(hi)) THEN
       v=0._dp; dv=0._dp; d2v=0._dp; RETURN
    ENDIF
    IF(x<r(1)) THEN
       v=f(1)+df(1)*(x-r(1)); dv=df(1); d2v=0._dp; RETURN
    ENDIF
    DO WHILE(hi-lo>1)
       mid=(lo+hi)/2
       IF(r(mid)>x) THEN
          hi=mid
       ELSE
          lo=mid
       ENDIF
    ENDDO
    h=r(hi)-r(lo); t=(x-r(lo))/h
    a=f(lo); b=h*df(lo)
    c=3._dp*(f(hi)-f(lo))-h*(2._dp*df(lo)+df(hi))
    d=2._dp*(f(lo)-f(hi))+h*(df(lo)+df(hi))
    v=a+t*(b+t*(c+t*d)); dv=(b+t*(2._dp*c+3._dp*t*d))/h
    d2v=(2._dp*c+6._dp*t*d)/h**2
  END SUBROUTINE

  REAL(DP) FUNCTION tb09_weight(r,rc) RESULT(w)
    REAL(DP), INTENT(IN) :: r,rc
    REAL(DP) :: t
    t=MAX(0._dp,MIN(1._dp,(r/rc-0.7_dp)/0.3_dp))
    w=1._dp-t**3*(10._dp+t*(-15._dp+6._dp*t))
  END FUNCTION

  SUBROUTINE tb09_partition(theta,wout,weights)
    REAL(DP), INTENT(IN) :: theta(:)
    REAL(DP), INTENT(OUT) :: wout,weights(:)
    REAL(DP) :: total
    wout=PRODUCT(1._dp-theta); total=SUM(theta)
    weights=0._dp
    IF(total>0._dp) weights=(1._dp-wout)*theta/total
  END SUBROUTINE

  SUBROUTINE tb09_nearest(x,dx)
    USE cell_base, ONLY : alat,at,bg
    REAL(DP), INTENT(IN) :: x(3)
    REAL(DP), INTENT(OUT) :: dx(3)
    REAL(DP) :: s(3),trial(3),best,bound(3)
    INTEGER :: a,b,c,lo(3),hi(3)
    s=MATMUL(TRANSPOSE(bg),x/alat); s=s-ANINT(s)
    dx=alat*MATMUL(at,s); best=SUM(dx**2)
    ! Every improving image obeys |s_i+n_i| <= |b_i|*sqrt(best)/alat.
    ! This finite bound also works for skew cells; a fixed 27-image search
    ! does not, in general, guarantee the nearest Cartesian image.
    bound=SQRT(SUM(bg**2,DIM=1))*SQRT(best)/alat+16._dp*EPSILON(best)
    lo=CEILING(-s-bound); hi=FLOOR(-s+bound)
    DO c=lo(3),hi(3)
       DO b=lo(2),hi(2)
          DO a=lo(1),hi(1)
             trial=alat*MATMUL(at,s+REAL([a,b,c],dp))
             IF(SUM(trial**2)<best) THEN
                best=SUM(trial**2); dx=trial
             ENDIF
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE

  SUBROUTINE tb09_interpolate(x,f)
    USE cell_base, ONLY : alat,bg
    USE fft_base, ONLY : dfftp
    REAL(DP), INTENT(IN) :: x(3)
    REAL(DP), INTENT(OUT) :: f(7)
    REAL(DP) :: s(3),t(3),w(4,3),z
    INTEGER :: q(3),ns(3),i,j,k,ii,jj,kk
    ns=[dfftp%nr1,dfftp%nr2,dfftp%nr3]
    s=MODULO(MATMUL(TRANSPOSE(bg),x/alat),1._dp)*REAL(ns,dp)
    q=FLOOR(s); t=s-REAL(q,dp)
    DO i=1,3
       z=t(i)
       w(:,i)=[-0.5_dp*z+z*z-0.5_dp*z**3, 1._dp-2.5_dp*z*z+1.5_dp*z**3, &
                    0.5_dp*z+2._dp*z*z-1.5_dp*z**3, -0.5_dp*z*z+0.5_dp*z**3]
    ENDDO
    f=0._dp
    DO k=1,4
       kk=MODULO(q(3)+k-2,ns(3))+1
       DO j=1,4
          jj=MODULO(q(2)+j-2,ns(2))+1
          DO i=1,4
             ii=MODULO(q(1)+i-2,ns(1))+1
             f=f+w(i,1)*w(j,2)*w(k,3)*tb09_pw_grid(ii,jj,kk,:)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE

  SUBROUTINE tb09_build_atom(ia,bec)
    USE ions_base, ONLY : nat,ityp
    USE uspp_param, ONLY : upf,nh,nhm
    USE uspp, ONLY : indv,nhtolm,lpx,lpl,ap
    USE atom, ONLY : g=>rgrid
    USE constants, ONLY : fpi
    INTEGER, INTENT(IN) :: ia
    REAL(DP), INTENT(IN) :: bec(nhm*(nhm+1)/2,nat,1)
    TYPE(paw_info) :: i
    REAL(DP), ALLOCATABLE :: ae(:,:,:),ps(:,:,:),rr(:),ur(:,:),tk(:,:)
    INTEGER :: m,l2,lm,side,nt
    REAL(DP) :: sgn
    i%a=ia; i%t=ityp(ia); i%m=g(i%t)%mesh; i%b=upf(i%t)%nbeta; i%l=upf(i%t)%lmax_rho+1
    m=i%m; l2=i%l**2; nt=i%t
    ASSOCIATE(c=>tb09_atoms(ia))
    c%radius=g(nt)%r(upf(nt)%paw%iraug)
    ALLOCATE(c%dn(m,l2),c%dnr(m,l2),c%dt(m,l2),c%dtr(m,l2))
    ALLOCATE(c%core(m),c%corer(m),c%tcore(m),c%tcorer(m))
    ALLOCATE(c%v_ae(m,l2),c%v_ps(m,l2),ae(m,l2,1),ps(m,l2,1),rr(m))
    rr=g(nt)%r
    CALL PAW_rho_lm(i,bec,upf(nt)%paw%pfunc,ae)
    CALL PAW_rho_lm(i,bec,upf(nt)%paw%ptfunc,ps,upf(nt)%qfuncl)
    c%dn=(ae(:,:,1)-ps(:,:,1))*SPREAD(g(nt)%rm2,2,l2)
    c%dn(:,1)=c%dn(:,1)+SQRT(fpi)*(upf(nt)%paw%ae_rho_atc-upf(nt)%rho_atc)
    c%dt=0._dp
    ALLOCATE(ur(m,i%b),tk(m,l2))
    DO side=1,2
       sgn=REAL(3-2*side,dp)
       IF(side==1) THEN
          ur=upf(nt)%aewfc*SPREAD(g(nt)%rm1,2,i%b)
       ELSE
          ur=upf(nt)%pswfc*SPREAD(g(nt)%rm1,2,i%b)
       ENDIF
       CALL tb09_tau_lm(rr,ur,bec(:,ia,1),indv(1:nh(nt),nt),nhtolm(1:nh(nt),nt), &
                        upf(nt)%lll,ap,lpx,lpl,tk)
       c%dt=c%dt+sgn*tk
    ENDDO
    c%dt(:,1)=c%dt(:,1)+SQRT(fpi)*(upf(nt)%paw%ae_tau_core-upf(nt)%paw%ps_tau_core)
    DO lm=1,l2
       CALL tb09_derivative(rr,c%dn(:,lm),c%dnr(:,lm))
       CALL tb09_derivative(rr,c%dt(:,lm),c%dtr(:,lm))
    ENDDO
    c%core=upf(nt)%rho_atc; c%tcore=upf(nt)%paw%ps_tau_core
    CALL tb09_derivative(rr,c%core,c%corer)
    CALL tb09_derivative(rr,c%tcore,c%tcorer)
    DEALLOCATE(ae,ps,rr,ur,tk)
    END ASSOCIATE
  END SUBROUTINE

  SUBROUTINE tb09_allfields(x,raw,ae,ps,wout,host_atom,host_weight)
    USE ions_base, ONLY : nat,ityp,tau
    USE cell_base, ONLY : alat,bg
    USE atom, ONLY : g=>rgrid
    INTEGER, OPTIONAL, INTENT(IN) :: host_atom
    REAL(DP), INTENT(IN) :: x(3),raw(6)
    REAL(DP), INTENT(OUT) :: ae(6),ps(6),wout
    REAL(DP), OPTIONAL, INTENT(OUT) :: host_weight
    REAL(DP) :: theta(nat),weights(nat),y0
    REAL(DP), ALLOCATABLE :: y(:,:),dy(:,:)
    REAL(DP) :: dx(3),r,r2(1),v,dv,d2v,t,dt,d2t,dr(1,3),frac(3),bound(3)
    INTEGER :: ia,nt,l2,lm,ip,ll,first,last
    ae=raw; ps=raw; wout=1._dp
    first=1; last=nat
    theta=0._dp
    ! For coreless atoms a known host quadrature needs only its overlapping
    ! neighbors. All field corrections are still summed before the nonlinear
    ! functional is evaluated. No isolated atomic ratios are added.
    DO ia=first,last
       nt=ityp(ia)
       IF(tb09_coreless.AND.PRESENT(host_atom)) THEN
          IF(.NOT.tb09_neighbors(ia,host_atom)) CYCLE
       ENDIF
       IF(tb09_coreless) THEN
          frac=MATMUL(TRANSPOSE(bg),x/alat-tau(:,ia)); frac=frac-ANINT(frac)
          bound=SQRT(SUM(bg**2,DIM=1))*tb09_atoms(ia)%radius/alat
          IF(ANY(ABS(frac)>bound+1.e-12_dp)) CYCLE
       ENDIF
       CALL tb09_nearest(x-alat*tau(:,ia),dx)
       r=MAX(SQRT(SUM(dx**2)),g(nt)%r(1))
       ASSOCIATE(c=>tb09_atoms(ia))
       IF(.NOT.tb09_coreless) THEN
       CALL tb09_radial(g(nt)%r,c%core,c%corer,r,v,dv,d2v)
       ps(1)=ps(1)+v; ps(2:4)=ps(2:4)+dv*dx/r; ps(5)=ps(5)+d2v+2._dp*dv/r
       CALL tb09_radial(g(nt)%r,c%tcore,c%tcorer,r,t,dt,d2t)
       ps(6)=ps(6)+t
       ENDIF
       IF(r<c%radius) THEN
          theta(ia)=tb09_weight(r,c%radius)
          l2=SIZE(c%dn,2)
          IF(l2==1) THEN
             y0=1._dp/SQRT(4._dp*ACOS(-1._dp))
             CALL tb09_radial(g(nt)%r,c%dn(:,1),c%dnr(:,1),r,v,dv,d2v)
             ae(1)=ae(1)+v*y0; ae(2:4)=ae(2:4)+dv*y0*dx/r
             ae(5)=ae(5)+(d2v+2._dp*dv/r)*y0
             CALL tb09_radial(g(nt)%r,c%dt(:,1),c%dtr(:,1),r,t,dt,d2t)
             ae(6)=ae(6)+t*y0
          ELSE
          ALLOCATE(y(1,l2),dy(1,l2)); dr(1,:)=dx; r2(1)=r*r
          CALL ylmr2(l2,1,dr,r2,y)
          DO lm=1,l2
             CALL tb09_radial(g(nt)%r,c%dn(:,lm),c%dnr(:,lm),r,v,dv,d2v)
             ae(1)=ae(1)+v*y(1,lm)
             ae(2:4)=ae(2:4)+dv*y(1,lm)*dx/r
             ll=INT(SQRT(REAL(lm-1,dp)))
             ae(5)=ae(5)+(d2v+2._dp*dv/r-REAL(ll*(ll+1),dp)*v/r**2)*y(1,lm)
             CALL tb09_radial(g(nt)%r,c%dt(:,lm),c%dtr(:,lm),r,t,dt,d2t)
             ae(6)=ae(6)+t*y(1,lm)
          ENDDO
          DO ip=1,3
             CALL dylmr2(l2,1,dr,r2,dy,ip)
             DO lm=1,l2
                CALL tb09_radial(g(nt)%r,c%dn(:,lm),c%dnr(:,lm),r,v,dv,d2v)
                ae(ip+1)=ae(ip+1)+v*dy(1,lm)
             ENDDO
          ENDDO
          DEALLOCATE(y,dy)
          ENDIF
       ENDIF
       END ASSOCIATE
    ENDDO
    ! Add the common pseudo-core contribution only once, after the deltas.
    ae=ae+(ps-raw)
    CALL tb09_partition(theta,wout,weights)
    IF(PRESENT(host_weight)) THEN
       IF(.NOT.PRESENT(host_atom)) CALL errore('tb09_allfields','Missing partition host',1)
       host_weight=weights(host_atom)
    ENDIF
  END SUBROUTINE

  SUBROUTINE tb09_model(fields,v)
    USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
    USE xc_lib, ONLY : xc_metagcx,xc
    USE constants, ONLY : e2
    REAL(DP), INTENT(IN) :: fields(:,:)
    REAL(DP), INTENT(OUT) :: v(:)
    REAL(DP), ALLOCATABLE :: n(:,:),gr(:,:,:),t(:,:),lap(:,:),ex(:),ec(:),v1(:,:),v2(:,:),v3(:,:), &
      c1(:,:),c2(:,:,:),c3(:,:),vx(:,:),vc(:,:)
    INTEGER :: m
    m=SIZE(fields,1)
    IF(ANY(.NOT.ieee_is_finite(fields))) CALL errore('tb09_model','Non-finite reconstructed field',1)
    IF(MINVAL(fields(:,1)) < -1.e-8_dp) &
       CALL errore('tb09_model','Negative reconstructed AE density; refine reconstruction',1)
    IF(MINVAL(fields(:,6)) < -1.e-8_dp) &
       CALL errore('tb09_model','Negative reconstructed AE kinetic density; refine reconstruction',1)
    ALLOCATE(n(m,1),gr(3,m,1),t(m,1),lap(m,1),ex(m),ec(m),v1(m,1),v2(m,1),v3(m,1), &
      c1(m,1),c2(1,m,1),c3(m,1),vx(m,1),vc(m,1))
    n(:,1)=MAX(fields(:,1),1.e-12_dp); gr(:,:,1)=TRANSPOSE(fields(:,2:4))
    t(:,1)=MAX(fields(:,6),1.e-14_dp); lap(:,1)=fields(:,5)
    CALL xc_metagcx(m,1,1,n,gr,t,ex,ec,v1,v2,v3,c1,c2,c3,lapl_=lap)
    v=e2*(v1(:,1)+c1(:,1))
    CALL xc(m,1,1,n,ex,ec,vx,vc)
    v=v+e2*(vx(:,1)+vc(:,1))
  END SUBROUTINE

  SUBROUTINE tb09_lda_reference(density,v)
    USE xc_lib, ONLY : xc
    USE constants, ONLY : e2,pi
    REAL(DP), INTENT(IN) :: density(:)
    REAL(DP), INTENT(OUT) :: v(:)
    REAL(DP) :: n(SIZE(density),1), ex(SIZE(density)),ec(SIZE(density)), &
                vx(SIZE(density),1),vc(SIZE(density),1)
    n(:,1)=MAX(ABS(density),1.e-8_dp)
    CALL xc(SIZE(density),1,1,n,ex,ec,vx,vc)
    v=e2*(-(3._dp*n(:,1)/pi)**(1._dp/3._dp)+vc(:,1))
  END SUBROUTINE

  SUBROUTINE PAW_tb09_prepare(rho,v,etxc,vtxc)
    USE scf, ONLY : scf_type
    USE ions_base, ONLY : nat,ityp,tau
    USE cell_base, ONLY : alat,at,bg,omega
    USE atom, ONLY : g=>rgrid
    USE uspp_param, ONLY : upf
    USE fft_base, ONLY : dfftp
    USE fft_types, ONLY : fft_index_to_3d
    USE scatter_mod, ONLY : gather_grid
    USE gvect, ONLY : gg,gvec=>g
    USE xc_lib, ONLY : set_libxc_ext_param
    USE mp, ONLY : mp_bcast,mp_sum
    USE io_global, ONLY : stdout
    USE constants, ONLY : e2
    TYPE(scf_type), INTENT(INOUT) :: rho
    REAL(DP), INTENT(OUT) :: v(:,:),etxc,vtxc
    REAL(DP), ALLOCATABLE :: local(:,:),full(:),gridfields(:,:),wg(:),vg(:),ref(:),psdensity(:),radfields(:,:), &
      vr(:),smooth(:),integ(:),volint(:),gl(:,:)
    REAL(DP) :: x(3),dx(3),raw(7),ae(6),ps(6),w,rr,stepvol,total,partvol,term,volterm,cpar,average
    INTEGER :: n1,n2,n3,ngrid,ix,iy,iz,ir,ia,ja,nt,m,nx,k,j,lm,q,np,ip,l2,image_bound(3)
    LOGICAL :: outside,overlap
    CALL start_clock('tb09_paw')
    IF(ALLOCATED(tb09_atoms)) DEALLOCATE(tb09_atoms)
    ALLOCATE(tb09_atoms(nat))
    DO ia=1,nat
       CALL tb09_build_atom(ia,rho%bec)
    ENDDO
    tb09_coreless=.TRUE.
    DO ia=1,nat
       IF(ANY(tb09_atoms(ia)%core/=0._dp).OR.ANY(tb09_atoms(ia)%tcore/=0._dp)) tb09_coreless=.FALSE.
    ENDDO
    IF(ALLOCATED(tb09_neighbors)) DEALLOCATE(tb09_neighbors)
    ALLOCATE(tb09_neighbors(nat,nat)); tb09_neighbors=.FALSE.; overlap=.FALSE.
    DO ia=1,nat
       tb09_neighbors(ia,ia)=.TRUE.
       DO ja=ia+1,nat
          CALL tb09_nearest(alat*(tau(:,ia)-tau(:,ja)),dx)
          IF(SQRT(SUM(dx**2))<tb09_atoms(ia)%radius+tb09_atoms(ja)%radius) THEN
             overlap=.TRUE.; tb09_neighbors(ia,ja)=.TRUE.; tb09_neighbors(ja,ia)=.TRUE.
             IF(.NOT.tb09_coreless.OR.ANY(ABS(upf(1:SIZE(upf))%zp-1._dp)>1.e-12_dp)) &
                CALL errore('PAW_tb09_prepare','Overlap experiment supports coreless hydrogen only',1)
             IF(SQRT(SUM(dx**2))<MAX(tb09_atoms(ia)%radius,tb09_atoms(ja)%radius)) &
                CALL errore('PAW_tb09_prepare','An augmentation sphere contains a neighboring nucleus',1)
          ENDIF
       ENDDO
    ENDDO
    IF(overlap) WRITE(stdout,'(5X,A)') &
       'TB09 PAW: experimental overlapping-sphere hydrogen reconstruction; accuracy must be checked.'
    DO ia=1,nat
       image_bound=FLOOR(2._dp*tb09_atoms(ia)%radius*SQRT(SUM(bg**2,DIM=1))/alat)
       DO iz=-image_bound(3),image_bound(3)
          DO iy=-image_bound(2),image_bound(2)
             DO ix=-image_bound(1),image_bound(1)
                IF(ix==0.AND.iy==0.AND.iz==0) CYCLE
                dx=alat*MATMUL(at,REAL([ix,iy,iz],dp))
                IF(SQRT(SUM(dx**2))<2._dp*tb09_atoms(ia)%radius) &
                   CALL errore('PAW_tb09_prepare','Self-overlapping augmentation sphere',1)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    n1=dfftp%nr1; n2=dfftp%nr2; n3=dfftp%nr3; ngrid=n1*n2*n3; stepvol=omega/REAL(ngrid,dp)
    IF(ALLOCATED(tb09_pw_grid)) DEALLOCATE(tb09_pw_grid)
    ALLOCATE(tb09_pw_grid(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x,7))
    ALLOCATE(local(dfftp%nnr,6),full(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x),gl(3,dfftp%nnr))
    local(:,1)=rho%of_r(:,1)
    CALL fft_gradient_g2r(dfftp,rho%of_g(:,1),gvec,gl)
    local(:,2:4)=TRANSPOSE(gl)
    CALL fft_laplacian(dfftp,rho%of_r(:,1),gg,local(:,5))
    local(:,6)=rho%kin_r(:,1)/e2
    tb09_pw_grid=0._dp
    DO ip=1,6
#if defined(__MPI)
       CALL gather_grid(dfftp,local(:,ip),full)
       CALL mp_bcast(full,dfftp%root,dfftp%comm)
#else
       full=local(:,ip)
#endif
       tb09_pw_grid(:,:,:,ip)=RESHAPE(full,SHAPE(tb09_pw_grid(:,:,:,ip)))
    ENDDO
    ALLOCATE(gridfields(ngrid,6),wg(ngrid),vg(ngrid),ref(ngrid),psdensity(ngrid))
    q=0; total=0._dp; partvol=0._dp
    DO iz=0,n3-1
       DO iy=0,n2-1
          DO ix=0,n1-1
             q=q+1; x=alat*MATMUL(at,[REAL(ix,dp)/n1,REAL(iy,dp)/n2,REAL(iz,dp)/n3])
             raw=tb09_pw_grid(ix+1,iy+1,iz+1,:)
             CALL tb09_allfields(x,raw(1:6),ae,ps,w)
             wg(q)=w; psdensity(q)=ps(1)
             IF(w>1.e-12_dp) THEN
                total=total+stepvol*w*SQRT(SUM(ae(2:4)**2))/MAX(ae(1),1.e-12_dp)
             ENDIF
             partvol=partvol+stepvol*w
             gridfields(q,:)=ae
             ! The reference potential is regular inside the atomic partition.
             IF(w<1.e-12_dp) THEN
                gridfields(q,:)=0._dp; gridfields(q,1)=MAX(ABS(ps(1)),1.e-8_dp)
                gridfields(q,6)=1._dp
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! Atomic contribution to the nonlinear integral, evaluated on the full
    ! reconstructed crystal density rather than isolated one-center densities.
    DO ia=1,nat
       nt=ityp(ia); m=g(nt)%mesh; nx=rad(nt)%nx
       ALLOCATE(integ(m),volint(m))
       integ=0._dp; volint=0._dp
       DO j=1,nx
          dx=[rad(nt)%sin_th(j)*rad(nt)%cos_phi(j),rad(nt)%sin_th(j)*rad(nt)%sin_phi(j),rad(nt)%cos_th(j)]
          DO k=1,upf(nt)%paw%iraug
             rr=g(nt)%r(k); w=tb09_weight(rr,tb09_atoms(ia)%radius)
             IF(w<1.e-14_dp) CYCLE
             x=alat*tau(:,ia)+rr*dx
             CALL tb09_interpolate(x,raw)
             CALL tb09_allfields(x,raw(1:6),ae,ps,w,host_atom=ia,host_weight=term)
             w=term*rad(nt)%ww(j)*rr**2
             integ(k)=integ(k)+w*SQRT(SUM(ae(2:4)**2))/MAX(ae(1),1.e-12_dp)
             volint(k)=volint(k)+w
          ENDDO
       ENDDO
       CALL simpson(m,integ,g(nt)%rab,term)
       CALL simpson(m,volint,g(nt)%rab,volterm)
       total=total+term; partvol=partvol+volterm
       DEALLOCATE(integ,volint)
    ENDDO
    average=total/omega; cpar=-0.012_dp+1.023_dp*SQRT(average)
    CALL set_libxc_ext_param(5,0,cpar)
    WRITE(stdout,'(5X,"TB09 PAW c = ",F12.8," (reconstructed average = ",ES14.6,")")') cpar,average
    WRITE(stdout,'(5X,"TB09 PAW partition volume relative error = ",ES12.4)') (partvol-omega)/omega
    IF(ABS(partvol/omega-1._dp)>0.005_dp) &
       CALL errore('PAW_tb09_prepare','Refine FFT/angular grids: inaccurate partition volume',1)
    CALL tb09_model(gridfields,vg)
    CALL tb09_lda_reference(psdensity,ref)
    vg=wg*vg+(1._dp-wg)*ref
    q=0
    DO iz=1,n3
       DO iy=1,n2
          DO ix=1,n1
             q=q+1; tb09_pw_grid(ix,iy,iz,7)=vg(q)
          ENDDO
       ENDDO
    ENDDO
    v=0._dp; vtxc=0._dp; etxc=0._dp
    DO ir=1,dfftp%nnr
       CALL fft_index_to_3d(ir,dfftp,ix,iy,iz,outside)
       IF(outside) CYCLE
       v(ir,1)=tb09_pw_grid(ix+1,iy+1,iz+1,7)
       vtxc=vtxc+v(ir,1)*rho%of_r(ir,1)*stepvol
    ENDDO
    CALL mp_sum(vtxc,dfftp%comm)
    ! AE and actual smooth-reference matrix elements of the local operator.
    DO ia=1,nat
       nt=ityp(ia); m=g(nt)%mesh; nx=rad(nt)%nx; l2=SIZE(tb09_atoms(ia)%dn,2)
       np=m*nx; ALLOCATE(radfields(np,6),vr(np),smooth(np)); radfields=0._dp; smooth=0._dp
       DO j=1,nx
          dx=[rad(nt)%sin_th(j)*rad(nt)%cos_phi(j),rad(nt)%sin_th(j)*rad(nt)%sin_phi(j),rad(nt)%cos_th(j)]
          DO k=1,upf(nt)%paw%iraug
             q=k+(j-1)*m; x=alat*tau(:,ia)+g(nt)%r(k)*dx
             CALL tb09_interpolate(x,raw)
             CALL tb09_allfields(x,raw(1:6),ae,ps,w,host_atom=ia)
             radfields(q,:)=ae; smooth(q)=raw(7)
          ENDDO
       ENDDO
       CALL tb09_model(radfields,vr)
       tb09_atoms(ia)%v_ae=0._dp; tb09_atoms(ia)%v_ps=0._dp
       DO lm=1,l2
          DO j=1,nx
             DO k=1,upf(nt)%paw%iraug
                q=k+(j-1)*m
                tb09_atoms(ia)%v_ae(k,lm)=tb09_atoms(ia)%v_ae(k,lm)+rad(nt)%wwylm(j,lm)*vr(q)
                tb09_atoms(ia)%v_ps(k,lm)=tb09_atoms(ia)%v_ps(k,lm)+rad(nt)%wwylm(j,lm)*smooth(q)
             ENDDO
          ENDDO
       ENDDO
       DEALLOCATE(radfields,vr,smooth)
    ENDDO
    DEALLOCATE(local,full,gl,gridfields,wg,vg,ref,psdensity)
    CALL stop_clock('tb09_paw')
  END SUBROUTINE PAW_tb09_prepare
