c
c This code assumes the covariances C_i,j are in an ASCII file called
c    cov
c with the first row listing the 11 covariances C_1,j, j=1,...,11
c the second row listing the 10 covariances C_2,j, j=2,...,11
c etc until the 11th row containing just the single scalar C_11,11.
c
c The first output (ASCII) file is called
c    eigenvalues
c and consists of two columns, the first listing the eigenvalues
c themselves (i.e. the variances of the PCs), NOT NECESSARILY IN
c DECREASING ORDER, and the second listing their proportional
c contribution to the total variance in % (and because the PCs
c are uncorrelated, the total variance is the sum of the individual
c variances).
c
c The second output (ASCII) file is called 
c       Wmatrix 
c and it contains the entries of the matrix W such that 
c              covariance = W * diagonal * Wtranspose
c so that the covariance of Wtranspose X is diagonal (if
c X was the original vector random variable),
c i.e. column 1 of W is the coefficients of the 1st PC,
c      column 2 of W is the coefficients of the 2nd PC,
c      column 3 of W is the coefficients of the 3rd PC,
c      column 4 of W is the coefficients of the 4th PC,
c      ...
c      column 11 of W is the coefficients of the 11th PC.
c The format of the Wmatrix file is 11 columns of 11 rows,
c showing  Wij  in row i , column j
c
c This code also writes to the screen the first three
c PCs

      double precision u(11,11), inv0(11,11)
      double precision inv(11,11),v(11,11),w(11)
c u(11,11) is the W matrix, and v(11,11) should turn out to be its transpose
c w(11) is the array of eigenvalues

      open(9,file='cov',status='old')
      read(9,*)inv0(1,1),inv0(1,2),inv0(1,3),inv0(1,4),inv0(1,5),
     . inv0(1,6),inv0(1,7),inv0(1,8),inv0(1,9),inv0(1,10),inv0(1,11)
      read(9,*)inv0(2,2),inv0(2,3),inv0(2,4),inv0(2,5),inv0(2,6),
     . inv0(2,7),inv0(2,8),inv0(2,9),inv0(2,10),inv0(2,11)
      read(9,*)inv0(3,3),inv0(3,4),inv0(3,5),inv0(3,6),inv0(3,7),
     . inv0(3,8),inv0(3,9),inv0(3,10),inv0(3,11)
      read(9,*)inv0(4,4),inv0(4,5),inv0(4,6),inv0(4,7),inv0(4,8),
     . inv0(4,9),inv0(4,10),inv0(4,11)
      read(9,*)inv0(5,5),inv0(5,6),inv0(5,7),inv0(5,8),inv0(5,9),
     . inv0(5,10),inv0(5,11)
      read(9,*)inv0(6,6),inv0(6,7),inv0(6,8),inv0(6,9),inv0(6,10),
     . inv0(6,11)
      read(9,*)inv0(7,7),inv0(7,8),inv0(7,9),inv0(7,10),inv0(7,11)
      read(9,*)inv0(8,8),inv0(8,9),inv0(8,10),inv0(8,11)
      read(9,*)inv0(9,9),inv0(9,10),inv0(9,11)
      read(9,*)inv0(10,10),inv0(10,11)
      read(9,*)inv0(11,11)
      close(9)

       do i=2,11
       do j=1,i-1
         inv0(i,j)=inv0(j,i)
       end do
       end do
c Put in inv0 the matrix you want to diagonalize. Then:

       do i=1,11
       do j=1,11
         u(i,j)=inv0(i,j)
       end do
       end do

      m=11
      n=11
      mp=11
      np=11
      call dsvdcmp(u,m,n,mp,np,w,v)

c So now inv0=u*w*(v-transp)
c and u should be equal to v if the original matrix was symmetric

      sum=0.
      do k=1,11
        sum=sum+w(k)
      end do
      open(9,file='eigenvalues',status='unknown')
c     write(9,*)'Eigenvalues:'
      do k=1,11
        write(9,*)w(k),w(k)*100./sum
      end do
      close(9)

      open(8,file='Wmatrix',status='unknown')
      do i=1,11
        write(8,*)(u(i,j),j=1,11)
        write(6,*)u(i,1),u(i,2),u(i,3)
      end do
      close(8)

      stop
      end

      function dpythag(a,b)
      double precision a,b,dpythag
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
      end

      subroutine dsvdcmp(a,m,n,mp,np,w,v)
      integer m,mp,n,np
      double precision a(11,11),v(11,11),w(11)
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
            h=dpythag(f,g)
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
          g=dpythag(f,1.0d0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0d0
          s=1.0d0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=dpythag(f,h)
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
            z=dpythag(f,h)
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
      end
