MODULE f_eof11

CONTAINS
SUBROUTINE eof11(a2cov, a1egval, a1cont, a2egvct)
!--------------------------------------------------------------------
! This code assumes the covariances C_i,j are in an ASCII file called
!    cov
! with the first row listing the 11 covariances C_1,j, j=1,...,11
! the second row listing the 10 covariances C_2,j, j=2,...,11
! etc until the 11th row containing just the single scalar C_11,11.
!
! The first output (ASCII) file is called
!    eigenvalues
! and consists of two columns, the first listing the eigenvalues
! themselves (i.e. the variances of the PCs), NOT NECESSARILY IN
! DECREASING ORDER, and the second listing their proportional
! contribution to the total variance in % (and because the PCs
! are uncorrelated, the total variance is the sum of the individual
! variances).
!
! The second output (ASCII) file is called 
!       Wmatrix 
! and it contains the entries of the matrix W such that 
!              covariance = W * diagonal * Wtranspose
! so that the covariance of Wtranspose X is diagonal (if
! X was the original vector random variable),
! i.e. column 1 of W is the coefficients of the 1st PC,
!      column 2 of W is the coefficients of the 2nd PC,
!      column 3 of W is the coefficients of the 3rd PC,
!      column 4 of W is the coefficients of the 4th PC,
!      ...
!      column 11 of W is the coefficients of the 11th PC.
! The format of the Wmatrix file is 11 columns of 11 rows,
! showing  Wij  in row i , column j
!
! This code also writes to the screen the first three
! PCs
!------------------------------------------------------------------
implicit none
!--- in -----------
double precision   a2cov(11,11)
!f2py intent(in)   a2cov

!--- out ----------
double precision   a1egval(11), a1cont(11)
!f2py intent(out)  a1egval, a1cont
double precision   a2egvct(11,11)
!f2py intent(out)  a2egvct

!--- calc  --------
double precision u(11,11), inv0(11,11)
double precision v(11,11),w(11)
double precision wsum
integer          m, n, mp, np
integer          i

! u(11,11) is the W matrix, and v(11,11) should turn out to be its transpose
! w(11) is the array of eigenvalues

       inv0 = a2cov

       !do i=2,11
       !do j=1,i-1
       !  inv0(i,j)=inv0(j,i)
       !end do
       !end do
! Put in inv0 the matrix you want to diagonalize. Then:

       u = inv0
       !do i=1,11
       !do j=1,11
       !  u(i,j)=inv0(i,j)
       !end do
       !end do

      m=11
      n=11
      mp=11
      np=11
      call dsvdcmp(u,m,n,mp,np,w,v)

! So now inv0=u*w*(v-transp)
! and u should be equal to v if the original matrix was symmetric

      wsum = sum(w)
      !wsum=0.
      !do k=1,11
      !  wsum=wsum+w(k)
      !end do
      !open(9,file='eigenvalues',status='unknown')
!     write(9,*)'Eigenvalues:'
      !do k=1,11
      !  write(9,*)w(k),w(k)*100./sum
      !end do
      !close(9)
      a1egval = w
      a1cont  = w*100./wsum


      !open(8,file='Wmatrix',status='unknown')
      do i=1,11
      !  write(8,*)(u(i,j),j=1,11)
      !  write(6,*)u(i,1),u(i,2),u(i,3)
      end do
      !close(8)
      a2egvct = u

      return
END SUBROUTINE eof11
!-----------------------------------------
function dpythag2(a,b)
implicit none
double precision  a,b
!f2py intent(in)  a,b
double precision  dpythag2
!f2py intent(out) dpythag2

double precision absa,absb
absa=abs(a)
absb=abs(b)
if(absa.gt.absb)then
  dpythag2=absa*sqrt(1.0d0+(absb/absa)**2)
else
  if(absb.eq.0.0d0)then
    dpythag2=0.0d0
  else
    dpythag2=absb*sqrt(1.0d0+(absa/absb)**2)
  endif
endif
return
end function dpythag2



function dpythag(a,b)
implicit none
double precision  a,b
!f2py intent(in)  a,b
double precision  dpythag
!f2py intent(out) dpythag

double precision absa,absb
absa=abs(a)
absb=abs(b)
if(absa.gt.absb)then
  dpythag=absa*sqrt(1.0d0+(absb/absa)**2)
else
  if(absb.eq.0.0d0)then
    dpythag=0.0d0
  else
    dpythag=absb*sqrt(1.0d0+(absa/absb)**2)
  endif
endif
return
end function dpythag

subroutine dsvdcmp(a,m,n,mp,np,w,v)
integer m,mp,n,np
!f2py intent(inout) m,mp,n,np
double precision a(11,11)
!f2py intent(inout) a
double precision v(11,11),w(11)
!f2py intent(inout) v, w


integer i,its,j,jj,k,l,nm
double precision anorm,c,f,g,h,s,scale,x,y,z,rv1(11),dpythag
      g=0.0d0
      scale=0.0d0
      anorm=0.0d0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if(i.le.m)then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
 11       continue
          if(scale.ne.0.0d0)then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
 12         continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do 15 j=l,n
              s=0.0d0
              do 13 k=i,m
                s=s+a(k,i)*a(k,j)
 13           continue
              f=s/h
              do 14 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
 14           continue
 15         continue
            do 16 k=i,m
              a(k,i)=scale*a(k,i)
 16         continue
          endif
        endif
        w(i)=scale *g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
 17       continue
          if(scale.ne.0.0d0)then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
 18         continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
 19         continue
            do 23 j=l,m
              s=0.0d0
              do 21 k=l,n
                s=s+a(j,k)*a(i,k)
 21           continue
              do 22 k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
 22           continue
 23         continue
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
 24         continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
 25   continue
      do 32 i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0d0)then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
 26         continue
            do 29 j=l,n
              s=0.0d0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
 27           continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
 28           continue
 29         continue
          endif
          do 31 j=l,n
            v(i,j)=0.0d0
            v(j,i)=0.0d0
 31       continue
        endif
        v(i,i)=1.0d0
        g=rv1(i)
        l=i
 32   continue
      do 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          a(i,j)=0.0d0
 33     continue
        if(g.ne.0.0d0)then
          g=1.0d0/g
          do 36 j=l,n
            s=0.0d0
            do 34 k=l,m
              s=s+a(k,i)*a(k,j)
 34         continue
            f=(s/a(i,i))*g
            do 35 k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
 35         continue
 36       continue
          do 37 j=i,m
            a(j,i)=a(j,i)*g
 37       continue
        else
          do 38 j= i,m
            a(j,i)=0.0d0
 38       continue
        endif
        a(i,i)=a(i,i)+1.0d0
 39   continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
 41       continue
 1        c=0.0d0
          s=1.0d0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)


            !print *,dpythag2(1.0d0,1.0d0)
            !print *,dpythag2(f,g)


            h=dpythag2(f,g)
            w(i)=h
            h=1.0d0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
 42         continue
 43       continue
 2        z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0d0)then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
 44           continue
            endif
            goto 3
          endif
          if(its.eq.30) write(6,*)'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
          g=dpythag2(f,1.0d0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0d0
          s=1.0d0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=dpythag2(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
 45         continue
            z=dpythag2(f,h)
            w(j)=z
            if(z.ne.0.0d0)then
              z=1.0d0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
 46         continue
 47       continue
          rv1(l)=0.0d0
          rv1(k)=f
          w(k)=x
 48     continue
 3      continue
 49   continue
return
end subroutine dsvdcmp

END MODULE f_eof11
