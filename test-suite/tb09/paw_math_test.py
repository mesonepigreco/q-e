#!/usr/bin/env python3
"""Compile the actual PAW radial helpers and test independent analytic fields.

These tests check a logarithmic mesh, a hydrogenic cusp, polynomial exactness,
and cubic Hermite derivatives. They do not use a second copy of the algorithm.
"""
from pathlib import Path
import argparse,subprocess,tempfile
import numpy as np
from numpy.polynomial.legendre import leggauss

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fc',default='gfortran');a=p.parse_args()
 root=Path(__file__).resolve().parents[2]
 source=(root/'PW/src/tb09_paw_reconstruction.f90').read_text()
 helpers=source[source.index('  SUBROUTINE tb09_derivative'):source.index('  SUBROUTINE tb09_nearest')]
 test='''
program fixtures
 use helpers_module
 implicit none
 integer,parameter :: m=1001
 real(dp) :: r(m),f(m),df(m),want(m),x,v,dv,d2v,err,err2
 integer :: i,k
 do i=1,m
   r(i)=1.e-4_dp*exp(real(i-1,dp)*log(2.e5_dp)/real(m-1,dp))
 enddo
 do k=0,4
   f=r**k
   want=0._dp
   if(k>0) want=k*r**(k-1)
   call tb09_derivative(r,f,df)
   err=maxval(abs(df-want)/max(1._dp,abs(want)))
   if(err>2.e-10_dp) stop 11
 enddo
 ! Independent hydrogenic n=exp(-2r): grad n=-2n, lap n=(4-4/r)n.
 f=exp(-2._dp*r)
 call tb09_derivative(r,f,df)
 err=0._dp;err2=0._dp
 do i=1,m-1
   x=(r(i)+r(i+1))/2._dp
   if(x<.02_dp.or.x>4._dp) cycle
   call tb09_radial(r,f,df,x,v,dv,d2v)
   err=max(err,abs(dv+2._dp*exp(-2._dp*x)))
   err2=max(err2,abs(d2v+2._dp*dv/x-(4._dp-4._dp/x)*exp(-2._dp*x)))
 enddo
 if(err>2.e-6_dp.or.err2>2.e-4_dp) then
   print *,err,err2
   stop 12
 endif
 if(abs(tb09_weight(0._dp,2._dp)-1._dp)>1.e-14_dp) stop 13
 if(abs(tb09_weight(2._dp,2._dp))>1.e-14_dp) stop 14
 print '(a,2es14.5)', 'analytic PAW radial fixtures PASS: ',err,err2
end program
'''
 # Independent Cartesian orbital: psi=exp(-r)*(a+b.x). Its gradient is
 # exp(-r)*(b-(a+b.x)*rhat), including all s-p cross terms. Obtain Gaunt
 # products by polynomial angular quadrature, independently of QE ap tables.
 z,wz=leggauss(8);phi=np.arange(16)*2*np.pi/16
 dirs=np.array([[np.sqrt(1-t*t)*np.cos(p),np.sqrt(1-t*t)*np.sin(p),t] for t in z for p in phi])
 weights=np.repeat(wz,16)*2*np.pi/16
 x,y,z=dirs.T;pi=np.pi
 Y=np.array([np.ones(len(z))/np.sqrt(4*pi),np.sqrt(3/(4*pi))*z,-np.sqrt(3/(4*pi))*x,-np.sqrt(3/(4*pi))*y,
 np.sqrt(5/(16*pi))*(3*z*z-1),-np.sqrt(15/(4*pi))*x*z,-np.sqrt(15/(4*pi))*y*z,
 np.sqrt(15/(16*pi))*(x*x-y*y),np.sqrt(15/(4*pi))*x*y]).T
 gaunt=np.einsum('ql,qi,qj,q->lij',Y,Y[:,:4],Y[:,:4],weights)
 angular_test='''
program angular_fixture
 use helpers_module
 implicit none
 integer,parameter :: m=501,nq=128
 real(dp) :: r(m),ur(m,2),tk(m,9),gaunt(9,4,4),bec(10),coeff(4),dir(nq,3),ylm(nq,9)
 real(dp) :: b(3),grad(3),actual,want,err,pi
 integer :: i,j,k,q,ij,lpx(4,4),lpl(4,4,9),indv(4),nhtolm(4),lll(2)
 pi=acos(-1._dp)
 open(10,file='gaunt.dat');read(10,*) gaunt;close(10)
 open(10,file='dirs.dat');read(10,*) dir;close(10)
 open(10,file='ylm.dat');read(10,*) ylm;close(10)
 coeff=[.7_dp,.2_dp,-.4_dp,.3_dp];b=[-coeff(3),-coeff(4),coeff(2)]
 ij=0
 do i=1,4
   do j=i,4
     ij=ij+1;bec(ij)=coeff(i)*coeff(j)
     if(i/=j) bec(ij)=2._dp*bec(ij)
   enddo
 enddo
 lpx=9
 do i=1,9
   lpl(:,:,i)=i
 enddo
 indv=[1,2,2,2];nhtolm=[1,2,3,4];lll=[0,1]
 do i=1,m
   r(i)=.001_dp*exp(real(i-1,dp)*log(4000._dp)/real(m-1,dp))
 enddo
 ur(:,1)=sqrt(4._dp*pi)*exp(-r)
 ur(:,2)=sqrt(4._dp*pi/3._dp)*r*exp(-r)
 call tb09_tau_lm(r,ur,bec,indv,nhtolm,lll,gaunt,lpx,lpl,tk)
 err=0._dp
 do q=1,nq
   do k=3,m-2
     grad=exp(-r(k))*(b-(coeff(1)+r(k)*sum(b*dir(q,:)))*dir(q,:))
     want=.5_dp*sum(grad**2)
     actual=sum(tk(k,:)*ylm(q,:))
     err=max(err,abs(actual-want))
   enddo
 enddo
 if(err>2.e-8_dp) then
   print *,err
   stop 21
 endif
 print '(a,es14.5)', 'Cartesian s-p kinetic-density fixture PASS: ',err
end program
'''
 with tempfile.TemporaryDirectory(prefix='tb09-paw-math-') as d:
  d=Path(d);src=d/'fixtures.f90';src.write_text('module helpers_module\nuse iso_fortran_env, only: dp=>real64\nimplicit none\ncontains\n'+helpers+'\nend module\n'+test)
  subprocess.run([a.fc,'-O2','-fcheck=all','-ffpe-trap=invalid,zero,overflow',str(src),'-o',str(d/'check')],cwd=d,check=True)
  subprocess.run([str(d/'check')],cwd=d,check=True)
  for name,data in [('gaunt',gaunt),('dirs',dirs),('ylm',Y)]:
   np.savetxt(d/(name+'.dat'),np.ravel(data,order='F'))
  src.write_text('module helpers_module\nuse iso_fortran_env, only: dp=>real64\nimplicit none\ncontains\n'+helpers+'\nend module\n'+angular_test)
  subprocess.run([a.fc,'-O2','-fcheck=all','-ffpe-trap=invalid,zero,overflow',str(src),'-o',str(d/'angular')],cwd=d,check=True)
  subprocess.run([str(d/'angular')],cwd=d,check=True)
  # Compare the actual nearest-image helper with exhaustive Cartesian search
  # on a strongly skew lattice, for which a fixed 27-image search can fail.
  matrix=np.array([[1.,3.2,1.9],[0.,1.,2.3],[0.,0.,1.]])
  np.savetxt(d/'at.dat',matrix.ravel(order='F'))
  np.savetxt(d/'bg.dat',np.linalg.inv(matrix).T.ravel(order='F'))
  nearest=source[source.index('  SUBROUTINE tb09_nearest'):source.index('  SUBROUTINE tb09_interpolate')]
  fixture='''program images
 use helpers_module
 use cell_base
 implicit none
 real(dp) :: frac(3),x(3),dx(3),trial(3),best,err
 integer :: q,i,j,k
 open(10,file='at.dat');read(10,*) at;close(10)
 open(10,file='bg.dat');read(10,*) bg;close(10)
 err=0._dp
 do q=1,100
   frac=.49_dp*sin([real(q,dp),real(2*q+3,dp),real(3*q+1,dp)])
   x=matmul(at,frac)
   call tb09_nearest(x,dx)
   best=huge(1._dp)
   do k=-12,12
     do j=-12,12
       do i=-12,12
         trial=x+matmul(at,real([i,j,k],dp));best=min(best,sum(trial**2))
       enddo
     enddo
   enddo
   err=max(err,abs(sum(dx**2)-best))
 enddo
 if(err>2.e-12_dp) stop 31
 print '(a,es14.5)', 'skew-cell exhaustive nearest-image fixture PASS: ',err
end program
'''
  src.write_text('module cell_base\nuse iso_fortran_env, only: dp=>real64\nreal(dp) :: alat=1._dp,at(3,3),bg(3,3)\nend module\nmodule helpers_module\nuse iso_fortran_env, only: dp=>real64\nimplicit none\ncontains\n'+nearest+'\nend module\n'+fixture)
  subprocess.run([a.fc,'-O2','-fcheck=all',str(src),'-o',str(d/'images')],cwd=d,check=True)
  subprocess.run([str(d/'images')],cwd=d,check=True)

if __name__=='__main__':main()
