! Fortran module for compiling the solvopt Solver for
! Local Nonlinear Optimization Problems used for Gubbin's PSD;
!  a colloid sim GSD file
! NOTE: requires compilation with compile-module
! NOTE: is run from a matching module_analysis fortran file
!       and sim-analysis Python script
! (Deepak Mangal and Rob Campbell)

subroutine solvopt(n,x,f,fun,flg,grad,flfc,func,flgc,gradc)
      
      implicit none

      character*(*) errmes,wrnmes,error2,error32,error42,error43
      character*(*) error52,error62,error63,error5,error6
      character*(*) warn1,warn20,warn21,warn4,warn31,warn32
      character*(*) warn09,warn08,termwarn0,termwarn1,appwarn
      character*(*) endwarn1, endwarn2, endwarn3, endwarn4
      parameter (errmes='SolvOpt error:')
      parameter (wrnmes='SolvOpt warning:')
      parameter (error2='Improper space dimension.')
      parameter (error32='Function equals infinity at the point.')
      parameter (error42='Gradient equals infinity at the starting point.',&
                &error43='Gradient equals zero at the starting point.')
      parameter (error52='<FUNC> returns infinite value at the point.',&
                &error62='<GRADC> returns infinite vector at the point.',&
                &error63='<GRADC>returns zero vector at an infeasible point.')
      parameter (error5='Function is unbounded.',&
                 &error6='Choose another starting point.')
      parameter (warn1='Gradient is zero, but stopping criteria are not fulfilled.',&
                 &warn20='Normal re-setting of a transformation matrix.',&
                 &warn21='Re-setting due to the use of a new penalty coefficient.')
      parameter (warn4='Iterations limit exceeded.',&
                 &warn31='The function is flat in certain directions.',&
                 &warn32='Trying to recover by shifting insensitive variables.',&
                 &warn09='Re-run from recorded point.',&
                 &warn08='Ravine with a flat bottom is detected.')
       parameter (termwarn0='SolvOpt: Normal termination.',&
                 &termwarn1='SolvOpt: Termination warning:',&
                 &appwarn='The above warning may be reasoned by inaccurate gradient approximation',&
                 &endwarn1='Premature stop is possible.Try to re-run the routine from the obtained point.',&
                 &endwarn2='Result may not provide the optimum. The function apparently has many extremum points.',&
                 &endwarn3='Result may be inaccurate in the coordinates. The function is flat at the solution.',&
                 &endwarn4='Stopping criteria are not fulfilled. The function is very steep at the solution.')
      
      logical:: flg,flgc,flfc, constr, app, appconstr
      logical:: FsbPnt, FsbPnt1, termflag, stopf
      logical:: stopping, dispwarn, Reset, ksm,knan,obj
      integer:: n, kstore, ajp,ajpp,knorms, k, kcheck, numelem
      integer:: dispdata, ld, mxtc, termx, limxterm, nzero, krerun
      integer:: warnno, kflat, stepvanish, i,j,ni,ii, kd,kj,kc,ip
      integer:: iterlimit, kg,k1,k2, kless,   allocerr
      real(8):: options(13),doptions(13) 
      real(8):: x(n),f
      real(8):: nsteps(3), gnorms(10), kk, nx
      real(8):: ajb,ajs, des, dq,du20,du10,du03
      real(8):: n_float, cnteps
      real(8):: low_bound, ZeroGrad, ddx, y
      real(8):: lowxbound, lowfbound, detfr, detxr, grbnd
      real(8):: fp,fp1,fc,f1,f2,fm,fopt,frec,fst, fp_rate
      real(8):: PenCoef, PenCoefNew
      real(8):: gamma,w,wdef,h1,h,hp
      real(8):: dx,ng,ngc,nng,ngt,nrmz,ng1,d,dd, laststep
      real(8):: zero,one,two,three,four,five,six,seven
      real(8):: eight,nine,ten,hundr
      real(8):: infty, epsnorm,epsnorm2,powerm12
      real(8),allocatable ::  B(:,:)
      real(8),allocatable ::  g(:)
      real(8),allocatable ::  g0(:)
      real(8),allocatable ::  g1(:)
      real(8),allocatable ::  gt(:)
      real(8),allocatable ::  gc(:)
      real(8),allocatable ::  z (:)
      real(8),allocatable ::  x1(:)
      real(8),allocatable ::  xopt(:)
      real(8),allocatable ::  xrec(:)
      real(8),allocatable ::  grec(:)
      real(8),allocatable ::  xx(:)
      real(8),allocatable ::  deltax(:)
      integer,allocatable ::   idx(:)
      character :: endwarn*100, allocerrstr*19
      external fun,grad,func,gradc
      
      data zero/0.d0/, one/1.d0/, two/2.d0/, three/3.d0/, four/4.d0/,five/5.d0/,&
           &six/6.d0/, seven/7.d0/, eight/8.d0/, nine/9.d0/,ten/1.d1/,&
           &hundr/1.d2/, powerm12/1.d-12/,infty /1.d100/, epsnorm /1.d-15/,&
           & epsnorm2 /1.d-30/,allocerrstr/'Allocation Error = '/

       options(1)=-1.d0
       options(2)=1.d-4
       options(3)=1.d-6
       options(4)=15.d3
       options(5)=-1.d0
       options(6)=1.d-8
       options(7)=2.5d0
       options(8)=1.d-11
       options(9)=0.d0
       options(10)=0.d0
       options(11)=0.d0
       options(12)=0.d0
       options(13)=0.d0

      if (n.lt.2) then
          print *, errmes
          print *, error2
        options(9)=-one
        goto 999
      endif  
      n_float=dble(n)

      allocate (B(n,n),stat=allocerr)
      if (allocerr.ne.0) then
         options(9)=-one
         print *,allocerrstr,allocerr
      endif   
      allocate (g(n),stat=allocerr)
      if (allocerr.ne.0) then
         options(9)=-one
         print *,allocerrstr,allocerr
      endif   
      allocate (g0(n),stat=allocerr)
      if (allocerr.ne.0) then
         options(9)=-one
         print *,allocerrstr,allocerr
      endif   
      allocate (g1(n),stat=allocerr)
      if (allocerr.ne.0) then
         options(9)=-one
         print *,allocerrstr,allocerr
      endif   
      allocate (gt(n),stat=allocerr)
      if (allocerr.ne.0) then
         options(9)=-one
         print *,allocerrstr,allocerr
      endif   
      allocate (gc(n),stat=allocerr)
      if (allocerr.ne.0) then
         options(9)=-one
         print *,allocerrstr,allocerr
      endif   
      allocate (z(n),stat=allocerr)
      if (allocerr.ne.0) then
         options(9)=-one
         print *,allocerrstr,allocerr
      endif   
      allocate (x1(n),stat=allocerr)
      if (allocerr.ne.0) then
         options(9)=-one
         print *,allocerrstr,allocerr
      endif   
      allocate (xopt(n),stat=allocerr)
      if (allocerr.ne.0) then
         options(9)=-one
         print *,allocerrstr,allocerr
      endif   
      allocate (xrec(n),stat=allocerr)
      if (allocerr.ne.0) then
         options(9)=-one
         print *,allocerrstr,allocerr
      endif   
      allocate (grec(n),stat=allocerr)
      if (allocerr.ne.0) then
         options(9)=-one
         print *,allocerrstr,allocerr
      endif   
      allocate (xx(n),stat=allocerr)
      if (allocerr.ne.0) then
         options(9)=-one
         print *,allocerrstr,allocerr
      endif   
      allocate (deltax(n),stat=allocerr)
      if (allocerr.ne.0) then
         options(9)=-one
         print *,allocerrstr,allocerr
      endif   
      allocate (idx(n),stat=allocerr)
      if (allocerr.ne.0) then
         options(9)=-one
         print *,allocerrstr,allocerr
      endif   


      app=.not.flg
      constr=flfc
      appconstr=.not.flgc

      call soptions(doptions)
      do i=1,8
            if (options(i).eq.zero) then
               options(i)=doptions(i)
            elseif (i.eq.2.or.i.eq.3.or.i.eq.6) then
               options(i)=dmax1(options(i),powerm12)
               options(i)=dmin1(options(i),one)
               if (i.eq.2)options(i)=dmax1(options(i),options(8)*hundr)
            elseif (i.eq.7) then
               options(7)=dmax1(options(i),1.5d0)
            endif
      enddo
               

                        
      options(10)=zero    !! counter for function calculations 
      options(11)=zero    !! counter for gradient calculations
      options(12)=zero    !! counter for constraint function calculations 
      options(13)=zero    !! counter for constraint gradient calculations
      iterlimit=idint(options(4))
      if (constr) then
        h1=-one           !! NLP: restricted to minimization 
        cnteps=options(6)
      else 
        h1=dsign(one,options(1))  !! Minimize resp. maximize a function
      endif
      k=0                         !! Iteration counter
      wdef=one/options(7)-one     !! Default space transf. coeff.


      ajb=one+1.d-1/n_float**2    !! Base I
      ajp=20  
      ajpp=ajp                    !! Start value for the power 
      ajs=1.15d0                  !! Base II
      knorms=0
      do i=1,10
       gnorms(i)=zero
      enddo  


      if (options(5).le.zero) then
         dispdata=0  
         if (options(5).eq.-one) then
            dispwarn=.false. 
         else 
            dispwarn=.true. 
         endif
      else 
         dispdata=idnint(options(5))
         dispwarn=.true.
      endif
      ld=dispdata



      dq=5.1d0           !! Step divider (at f_{i+1}>gamma*f_{i})
      du20=two
      du10=1.5d0
      du03=1.05d0        !! Step multipliers (at certain steps made)
      kstore=3
      do i=1,kstore
       nsteps(i)=zero    !! Steps made at the last 'kstore' iterations 
      enddo 
      if (app) then
        des=6.3d0        !! Desired number of steps per 1-D search
      else
        des=3.3d0
      endif
      mxtc=3             !! Number of trial cycles (steep wall detect)

      termx=0
      limxterm=50        !! Counter and limit for x-criterion
 
      ddx=dmax1(1.d-11,options(8))     

      low_bound=-one+1.d-4     !! Lower bound cosine used to detect aravine
      ZeroGrad=n_float*1.d-16  !! Lower bound for a gradient norm
      nzero=0                  !! Zero-gradient events counter

      lowxbound=dmax1(options(2),1.d-3)  

      lowfbound=options(3)**2
      krerun=0                 !! Re-run events counter
      detfr=options(3)*hundr   !! Relative error for f/f_{record}
      detxr=options(2)*ten     !! Relative error fornorm(x)/norm(x_{record})
      warnno=0                 !! the number of warn.mess. to end with
      kflat=0                  !! counter for points of flatness
      stepvanish=0             !! counter for vanished steps
      stopf=.false.



      call fun(x,f)
      options(10)=options(10)+one
      if (dabs(f).ge.infty) then
         if (dispwarn) then
            print *,errmes
            print *,error32
            print *,error6
         endif   
         options(9)=-three
         goto 999
      endif
      do i=1,n
        xrec(i)=x(i)
      enddo  
      frec=f     !! record point and function value
   
      if (constr)  then
          kless=0
          fp=f
          call func(x,fc)
          options(12)=options(12)+one 
          if (dabs(fc).ge.infty) then
             if (dispwarn) then
                print *,errmes
                print *,error52
                print *,error6
             endif
             options(9)=-five
             goto 999
          endif   
        PenCoef=one          !! first rough approximation
        if (fc.le.cnteps) then  
         FsbPnt=.true.       !! feasible point  
         fc=zero             
        else
         FsbPnt=.false. 
        endif
        f=f+PenCoef*fc
      endif   

      if (app) then
        do i=1,n
         deltax(i)=h1*ddx
        enddo
        obj=.true.
        if (constr) then
           call apprgrdn(n,g,x,fp,fun,deltax,obj)
        else
           call apprgrdn(n,g,x,f,fun,deltax,obj)  
        endif
        options(10)=options(10)+n_float
      else
        call grad(x,g)
        options(11)=options(11)+one
      endif
      ng=zero      
      do i=1,n
         ng=ng+g(i)*g(i)
      enddo
      ng=dsqrt(ng)
      if (ng.ge.infty) then
         if (dispwarn) then
            print *,errmes
            print *,error42
            print *,error6
         endif
         options(9)=-four
         goto 999
      elseif (ng.lt.ZeroGrad) then
         if (dispwarn) then
            print *,errmes
            print *,error43
            print *,error6
         endif
         options(9)=-four
         goto 999
      endif
      if (constr) then
       if (.not.FsbPnt) then
         if (appconstr) then
            do j=1,n
              if (x(j).ge.zero) then
                 deltax(j)=ddx
              else
                 deltax(j)=-ddx
              endif
            enddo
            obj=.false.     
            call apprgrdn(n,gc,x,fc,func,deltax,obj)
         else
            call gradc(x,gc)  
         endif
         ngc=zero      
         do i=1,n
           ngc=ngc+gc(i)*gc(i)
         enddo
         ngc=dsqrt(ngc)
         if (ng.ge.infty) then
            if (dispwarn) then
               print *,errmes
               print *,error62
               print *,error6
            endif
            options(9)=-six
            goto 999
         elseif (ng.lt.ZeroGrad) then
            if (dispwarn) then
               print *,errmes
               print *,error63
            endif
            options(9)=-six
            goto 999
         endif
         do i=1,n
           g(i)=g(i)+PenCoef*gc(i)
         enddo
         ng=zero      
         do i=1,n
           ng=ng+g(i)*g(i)
           grec(i)=g(i) 
         enddo
         ng=dsqrt(ng)
       endif
      endif
      do i=1,n
        grec(i)=g(i) 
      enddo
      nng=ng


      d=zero
      do i=1,n
        if (d.lt.dabs(x(i))) d=dabs(x(i))
      enddo  
      h=h1*dsqrt(options(2))*d                  !! smallest possiblestepsize
      if (dabs(options(1)).ne.one) then 
        h=h1*dmax1(dabs(options(1)),dabs(h))    !! user-suppliedstepsize
      else  
          h=h1*dmax1(one/dlog(ng+1.1d0),dabs(h)) !! calculated stepsize
      endif


      do while (.true.)
        kcheck=0                       !! Set checkpoint counter.
        kg=0                           !! stepsizes stored
        kj=0                           !! ravine jump counter
        do i=1,n
          do j=1,n
            B(i,j)=zero
          enddo
          B(i,i)=one                   !! re-set transf. matrix toidentity 
          g1(i)=g(i)
        enddo     
        fst=f
        dx=0 
    


   
        do while (.true.)
          k=k+1
          kcheck=kcheck+1
          laststep=dx

           gamma=one+dmax1(ajb**((ajp-kcheck)*n),two*options(3))
           gamma=dmin1 ( gamma,ajs**dmax1(one,dlog10(nng+one)) )
      
       ngt=zero
       ng1=zero
       dd=zero
       do i=1,n
         d=zero
         do j=1,n
            d=d+B(j,i)*g(j)
         enddo
         gt(i)=d
         dd=dd+d*g1(i)
         ngt=ngt+d*d
         ng1=ng1+g1(i)*g1(i)
       enddo
       ngt=dsqrt(ngt)      
       ng1=dsqrt(ng1)
       dd=dd/ngt/ng1      
       
       w=wdef       
      
       if (dd.lt.low_bound) then
        if (kj.eq.2) then
          do i=1,n
           xx(i)=x(i)
          enddo
        endif   
        if (kj.eq.0) kd=4
        kj=kj+1
        w=-.9d0              !! use large coef. of space dilation 
        h=h*two  
        if (kj.gt.2*kd) then
          kd=kd+1
          warnno=1
          endwarn=endwarn1  
          do i=1,n
            if (dabs(x(i)-xx(i)).lt.epsnorm*dabs(x(i))) then
             if (dispwarn)  then
                print *,wrnmes
                print *,warn08
             endif
            endif
          enddo
        endif  
       else
        kj=0 
       endif

      
       nrmz=zero
       do i=1,n
         z(i)=gt(i)-g1(i)
         nrmz=nrmz+z(i)*z(i)
       enddo  
       nrmz=dsqrt(nrmz)
       if (nrmz.gt.epsnorm*ngt) then 
        do i=1,n
         z(i)=z(i)/nrmz               
        enddo 


        d = zero
        do i=1,n
          d=d+z(i)*gt(i)
        enddo
        ng1=zero
        d = d*w
        do i=1,n
          dd=zero
          g1(i)=gt(i)+d*z(i)
          ng1=ng1+g1(i)*g1(i)
          do j=1,n
             dd=dd+B(i,j)*z(j)
          enddo 
          dd=w*dd
          do j=1,n
            B(i,j)=B(i,j)+dd*z(j)
          enddo
        enddo      
        ng1=dsqrt(ng1)
       else
        do i=1,n
         z(i)=zero
         g1(i)=gt(i)  
        enddo
        nrmz=zero
       endif
       do i=1,n
           gt(i)=g1(i)/ng1
       enddo
        do i=1,n
          d=zero
            do j=1,n
               d=d+B(i,j)*gt(j)
            enddo  
          g0(i)=d
        enddo


        if (kcheck.gt.1) then
           numelem=0
           do i=1,n
              if (dabs(g(i)).gt.ZeroGrad) then
                 numelem=numelem+1
                 idx(numelem)=i
              endif   
           enddo   
           if (numelem.gt.0) then
              grbnd=epsnorm*dble(numelem**2)
              ii=0
              do i=1,numelem
                 j=idx(i)
                 if (dabs(g1(j)).le.dabs(g(j))*grbnd) ii=ii+1
              enddo   
              if (ii.eq.n .or. nrmz.eq.zero) then
                if (dispwarn) then
                  print *,wrnmes
                  print *,warn20
                endif
                if (dabs(fst-f).lt.dabs(f)*1.d-2) then
                   ajp=ajp-10*n
                else
                   ajp=ajpp
                endif
                h=h1*dx/three
                k=k-1 
                exit
              endif
           endif 
        endif

 
        do i=1,n
         xopt(i)=x(i)
        enddo 
        fopt=f   
        k1=0
        k2=0
        ksm=.false.
        kc=0
        knan=.false.
        hp=h
        if (constr) Reset=.false.
 
        do while (.true.)
         do i=1,n
          x1(i)=x(i)
         enddo 
         f1=f   
         if (constr) then
           FsbPnt1=FsbPnt
           fp1=fp
         endif
         
         do i=1,n
            x(i)=x(i)+hp*g0(i)
         enddo
           ii=0
           do i=1,n   
            if (dabs(x(i)-x1(i)).lt.dabs(x(i))*epsnorm) ii=ii+1
           enddo    
         
         call fun(x,f)
         options(10)=options(10)+one 
         if (h1*f.ge.infty) then
            if (dispwarn) then
              print *,errmes 
              print *,error5
            endif
            options(9)=-seven
            goto 999
         endif
         if (constr) then
           fp=f
           call func(x,fc)
           options(12)=options(12)+one
           if (dabs(fc).ge.infty) then 
               if (dispwarn) then
                  print *,errmes
                  print *,error52
                  print *,error6
               endif
               options(9)=-five
               goto 999
           endif
           if (fc.le.cnteps) then
              FsbPnt=.true.
              fc=zero
           else
              FsbPnt=.false.
              fp_rate=fp-fp1 
              if (fp_rate.lt.-epsnorm) then
               if (.not.FsbPnt1) then
                d=zero
                do i=1,n
                  d=d+(x(i)-x1(i))**2
                enddo
                d=dsqrt(d)  
                PenCoefNew=-1.5d1*fp_rate/d
                if (PenCoefNew.gt.1.2d0*PenCoef) then
                  PenCoef=PenCoefNew
                  Reset=.true.
                  kless=0 
                  f=f+PenCoef*fc
                  exit
                endif
               endif 
              endif 
           endif
           f=f+PenCoef*fc
         endif

         if (dabs(f).ge.infty) then
             if (dispwarn) then
               print *,wrnmes
               print *,error32
             endif
             if (ksm.or.kc.ge.mxtc) then
                options(9)=-three
                goto 999
             else 
                k2=k2+1
                k1=0 
                hp=hp/dq 
                do i=1,n
                 x(i)=x1(i)
                enddo 
                f=f1 
                knan=.true. 
                if (constr) then
                  FsbPnt=FsbPnt1 
                  fp=fp1 
                endif
             endif

         elseif (ii.eq.n) then
                stepvanish=stepvanish+1
                if (stepvanish.ge.5) then
                    options(9)=-ten-four
                    if (dispwarn) then 
                       print *,termwarn1
                       print *,endwarn4 
                    endif
                    goto 999
                else 
                    do i=1,n
                     x(i)=x1(i)
                    enddo 
                    f=f1 
                    hp=hp*ten 
                    ksm=.true.
                    if (constr) then
                       FsbPnt=FsbPnt1 
                       fp=fp1 
                    endif
                endif

         elseif (h1*f.lt.h1*gamma**idint(dsign(one,f1))*f1) then
             if (ksm) exit  
             k2=k2+1
             k1=0 
             hp=hp/dq
             do i=1,n
              x(i)=x1(i)
             enddo 
             f=f1 
             if (constr) then
                FsbPnt=FsbPnt1 
                fp=fp1 
             endif
             if (kc.ge.mxtc) exit

         else   
             if (h1*f.le.h1*f1) exit

             k1=k1+1 
             if (k2.gt.0) kc=kc+1 
             k2=0
             if (k1.ge.20) then 
                 hp=du20*hp 
             elseif (k1.ge.10) then 
                 hp=du10*hp
             elseif (k1.ge.3) then
                 hp=du03*hp
             endif
         endif
        enddo


        dx=zero
        do i=1,n
           dx=dx+(xopt(i)-x(i))**2
        enddo
        dx=dsqrt(dx)   
        if (kg.lt.kstore)  kg=kg+1
        if (kg.ge.2)  then
           do i=kg,2,-1
             nsteps(i)=nsteps(i-1)
           enddo
        endif
        d=zero   
        do i=1,n
           d=d+g0(i)*g0(i)
        enddo
        d=dsqrt(d)
        nsteps(1)=dx/(dabs(h)*d)
        kk=zero
        d=zero
        do i=1,kg
           dd=dble(kg-i+1) 
           d=d+dd
           kk=kk+nsteps(i)*dd
        enddo
        kk=kk/d   
        if     (kk.gt.des) then
             if (kg.eq.1) then
                h=h*(kk-des+one)
             else   
                h=h*dsqrt(kk-des+one) 
             endif
        elseif (kk.lt.des) then
             h=h*dsqrt(kk/des)  
        endif

        if (ksm) stepvanish=stepvanish+1


        if (app) then
          do j=1,n
            if (g0(j).ge.zero) then
               deltax(j)=h1*ddx
            else
               deltax(j)=-h1*ddx
            endif
          enddo
          obj=.true.     
          if (constr)  then
             call apprgrdn(n,g,x,fp,fun,deltax,obj)
          else 
             call apprgrdn(n,g,x,f,fun,deltax,obj)
          endif
          options(10)=options(10)+n_float
        else
          call grad(x,g)
          options(11)=options(11)+one
        endif
        ng=zero
        do i=1,n
          ng=ng+g(i)*g(i)
        enddo
        ng=dsqrt(ng)
        if (ng.ge.infty) then
         if (dispwarn) then
           print *,errmes
           print *,error42
         endif
         options(9)=-four
         goto 999
        elseif (ng.lt.ZeroGrad) then
         if (dispwarn) then
           print *,wrnmes
           print *,warn1
         endif
         ng=ZeroGrad
        endif
      
        if (constr) then
         if (.not.FsbPnt) then
           if (ng.lt.1.d-2*PenCoef) then
              kless=kless+1
              if (kless.ge.20) then
                 PenCoef=PenCoef/ten
                 Reset=.true.
                 kless=0
              endif   
           else
              kless=0
           endif      
           if (appconstr) then
                 do j=1,n
                   if (x(j).ge.zero) then
                      deltax(j)=ddx
                   else
                      deltax(j)=-ddx
                   endif
                 enddo
                 obj=.false.     
                 call apprgrdn(n,gc,x,fc,func,deltax,obj)
                 options(12)=options(12)+n_float 
           else
                 call gradc(x,gc)  
                 options(13)=options(13)+one 
           endif
           ngc=zero
           do i=1,n
              ngc=ngc+gc(i)*gc(i)
           enddo
           ngc=dsqrt(ngc)
           if (ngc.ge.infty) then
                  if (dispwarn) then
                     print *,errmes
                     print *,error62
                  endif
                  options(9)=-six
                  goto 999
           elseif (ngc.lt.ZeroGrad .and. .not.appconstr) then 
                  if (dispwarn) then
                     print *,errmes
                     print *,error63
                  endif
                  options(9)=-six 
                  goto 999
           endif
           do i=1,n
             g(i)=g(i)+PenCoef*gc(i) 
           enddo
           ng=zero
           do i=1,n
              ng=ng+g(i)*g(i)
           enddo
           ng=dsqrt(ng)
           if (Reset) then
              if (dispwarn) then
                 print *,wrnmes
                 print *,warn21 
              endif
              h=h1*dx/three
              k=k-1
              nng=ng
              exit
           endif
         endif 
        endif
        if (h1*f.gt.h1*frec) then
          frec=f 
          do i=1,n
            xrec(i)=x(i)
            grec(i)=g(i)
          enddo 
        endif

       if (ng.gt.ZeroGrad) then
        if (knorms.lt.10)  knorms=knorms+1
        if (knorms.ge.2)  then
          do i=knorms,2,-1
           gnorms(i)=gnorms(i-1) 
          enddo
        endif  
        gnorms(1)=ng
        nng=one  
          do i=1,knorms
            nng=nng*gnorms(i)
          enddo
        nng=nng**(one/dble(knorms))
       endif

       nx=zero
       do i=1,n
        nx=nx+x(i)*x(i)
       enddo
       nx=dsqrt(nx)         
  

       if (k.eq.ld) then
         print *,'Iteration # ..... Function Value ..... ','Step Value ..... Gradient Norm'
         print '(5x,i5,7x,g13.5,6x,g13.5,7x,g13.5)', k,f,dx,ng
         ld=k+dispdata
       endif


      termflag=.true.
      if (constr) then
        if (.not.FsbPnt) termflag=.false.
      endif
      if(kcheck.le.5.or.kcheck.le.12.and.ng.gt.one)termflag=.false.
      if(kc.ge.mxtc .or. knan)termflag=.false.

       if (termflag) then
           ii=0
           stopping=.true.
           do i=1,n
             if (dabs(x(i)).ge.lowxbound) then
                ii=ii+1
                idx(ii)=i
                if (dabs(xopt(i)-x(i)).gt.options(2)*dabs(x(i))) then
                  stopping=.false.
                endif  
             endif
           enddo
           if (ii.eq.0 .or. stopping)  then
                stopping=.true.
                termx=termx+1
                d=zero
                do i=1,n
                  d=d+(x(i)-xrec(i))**2
                enddo
                d=dsqrt(d)  

                if(dabs(f-frec).gt.detfr*dabs(f) .and. dabs(f-fopt).le.options(3)*dabs(f) .and. krerun.le.3 .and. .not. constr) then
                   stopping=.false.
                   if (ii.gt.0) then
                    do i=1,ii
                     j=idx(i)
                     if (dabs(xrec(j)-x(j)).gt.detxr*dabs(x(j))) then
                       stopping=.true.
                       exit
                     endif
                    enddo
                   endif     
                   if (stopping) then
                      if (dispwarn) then
                        print *,wrnmes
                        print *,warn09
                      endif
                      ng=zero
                      do i=1,n
                       x(i)=xrec(i)
                       g(i)=grec(i)
                       ng=ng+g(i)*g(i)
                      enddo
                      ng=dsqrt(ng)  
                      f=frec
                      krerun=krerun+1
                      h=h1*dmax1(dx,detxr*nx)/dble(krerun)
                      warnno=2
                      endwarn=endwarn2
                      exit
                   else 
                      h=h*ten
                   endif
                elseif(dabs(f-frec).gt.options(3)*dabs(f) .and. d.lt.options(2)*nx .and. constr) then
                   continue
                elseif  (dabs(f-fopt).le.options(3)*dabs(f) .or. dabs(f).le.lowfbound .or. &
                        & (dabs(f-fopt).le.options(3).and. termx.ge.limxterm )) then
                  if (stopf) then
                   if (dx.le.laststep) then
                    if (warnno.eq.1 .and. ng.lt.dsqrt(options(3))) then
                       warnno=0 
                    endif
                    if (.not.app) then
                      do i=1,n
                       if (dabs(g(i)).le.epsnorm2) then
                         warnno=3
                         endwarn=endwarn3
                         exit
                       endif
                      enddo
                    endif  
                    if (warnno.ne.0) then
                       options(9)=-dble(warnno)-ten
                       if (dispwarn) then 
                         print *,termwarn1
                         print *,endwarn
                         if (app) print *,appwarn
                       endif    
                    else 
                       options(9)=dble(k)
                       if (dispwarn) print *,termwarn0
                    endif
                    goto 999
                   endif
                  else
                   stopf=.true.
                  endif 
                elseif (dx.lt.powerm12*dmax1(nx,one) .and. termx.ge.limxterm ) then
                     options(9)=-four-ten
                     if (dispwarn) then 
                       print *,termwarn1 
                       print *,endwarn4
                       if (app) print *,appwarn
                       f=frec
                       do i=1,n
                        x(i)=xrec(i)
                       enddo
                     endif
                     goto 999
                endif
           endif
       endif 

            if(k.eq.iterlimit) then
                options(9)=-nine
                if (dispwarn) then
                  print *,wrnmes
                  print *,warn4 
                endif
                goto 999
            endif


          if (constr) then 
            if (ng.le.ZeroGrad) then
                if (dispwarn) then  
                  print *,termwarn1
                  print *,warn1
                endif
                options(9)=-eight
                goto 999
            endif
          else  
            if (ng.le.ZeroGrad) then
             nzero=nzero+1
             if (dispwarn) then
               print *,wrnmes
               print *,warn1
             endif
             if (nzero.ge.3) then
               options(9)=-eight
               goto 999
             endif  
             do i=1,n
               g0(i)=-h*g0(i)/two
             enddo  
             do i=1,10
               do j=1,n
                x(j)=x(j)+g0(j)
               enddo                
               call fun(x,f)
               options(10)=options(10)+one
               if (dabs(f).ge.infty) then 
                 if (dispwarn) then
                   print *,errmes
                   print *,error32
                 endif
                 options(9)=-three
                 goto 999
               endif
               if (app) then
                   do j=1,n
                     if (g0(j).ge.zero) then
                        deltax(j)=h1*ddx
                     else
                        deltax(j)=-h1*ddx
                     endif
                   enddo
                   obj=.true.     
                   call apprgrdn(n,g,x,f,fun,deltax,obj)
                   options(10)=options(10)+n_float
               else
                   call grad(x,g)
                   options(11)=options(11)+one
               endif    
               ng=zero
               do j=1,n
                  ng=ng+g(j)*g(j)
               enddo
               ng=dsqrt(ng)
               if (ng.ge.infty) then
                    if (dispwarn) then
                      print *,errmes
                      print *,error42
                    endif
                    options(9)=-four
                    goto 999
               endif
               if (ng.gt.ZeroGrad) exit
             enddo
             if (ng.le.ZeroGrad) then
                if (dispwarn) then  
                  print *,termwarn1
                  print *,warn1
                endif
                options(9)=-eight
                goto 999
             endif
             h=h1*dx
             exit
            endif
          endif  


          if (.not.constr .and. dabs(f-fopt).lt.dabs(fopt)*options(3) .and. kcheck.gt.5  .and. ng.lt.one ) then
          
           ni=0
           do i=1,n
             if (dabs(g(i)).le.epsnorm2) then
               ni=ni+1
               idx(ni)=i
             endif
           enddo    
           if (ni.ge.1 .and. ni.le.n/2 .and. kflat.le.3) then
             kflat=kflat+1
             if (dispwarn) then
                print *,wrnmes
                print *,warn31
             endif
             warnno=1
             endwarn=endwarn1
             do i=1,n
               x1(i)=x(i)
             enddo  
             fm=f 
             do i=1,ni
              j=idx(i)
              f2=fm
              y=x(j)
              if (y.eq.zero) then
                x1(j)=one
              elseif (dabs(y).lt.one) then
                x1(j)=dsign(one,y)
              else
                x1(j)=y
              endif
              do ip=1,20
               x1(j)=x1(j)/1.15d0
               call fun(x1,f1)
               options(10)=options(10)+one
               if (dabs(f1).lt.infty) then
                 if (h1*f1.gt.h1*fm) then
                   y=x1(j)
                   fm=f1
                 elseif (h1*f2.gt.h1*f1) then
                   exit
                 elseif (f2.eq.f1) then
                   x1(j)=x1(j)/1.5d0   
                 endif  
                 f2=f1   
               endif
              enddo
              x1(j)=y
             enddo
             if (h1*fm.gt.h1*f) then
              if (app) then
                do j=1,n
                  deltax(j)=h1*ddx
                enddo
                obj=.true.     
                call apprgrdn(n,gt,x1,fm,fun,deltax,obj)
                options(10)=options(10)+n_float
              else
                call grad(x1,gt)
                options(11)=options(11)+one
              endif
              ngt=zero
              do i=1,n
                ngt=ngt+gt(i)*gt(i)
              enddo  
              if (ngt.gt.epsnorm2 .and. ngt.lt.infty) then  
                if (dispwarn) print *,warn32
                do i=1,n
                 x(i)=x1(i)
                 g(i)=gt(i)
                enddo  
                ng=ngt
                f=fm
                h=h1*dx/three
                options(3)=options(3)/five
                exit
              endif   !! regular gradient
             endif   !! a better value has been found
           endif   !! function is flat
          endif   !! pre-conditions are fulfilled

       enddo   !! iterations
      enddo   !! restart

999   continue
      

       deallocate (idx,deltax,xx,grec,xrec,xopt,x1,z,gc,gt,g1,g0,g,B)
      
      end 
      
      subroutine null
      end




      subroutine apprgrdn(n,g,x,f,fun,deltax,obj)

      implicit none
      double precision x(*), g(*), f, deltax(*)
      double precision lowbndobj,lowbndcnt,d,y,fi,one,ten,half
      integer n, i, j
      logical center,obj
      external fun
      data lowbndobj /2.d-10/, lowbndcnt /5.d-15/
      data one /1.d0/, ten /1.d1/, half /.5d0/
      do i=1,n
         y=x(i)
         d=dmax1(lowbndcnt,dabs(y))
         d=deltax(i)*d
         if (obj) then
           if (dabs(d).lt.lowbndobj) then
              d=lowbndobj*dsign(one,deltax(i))
              center=.true.
           else
              center=.false.
           endif
         else
           if (dabs(d).lt.lowbndcnt) then
              d=lowbndcnt*dsign(one,deltax(i))
           endif
         endif
         x(i)=y+d
         call fun(x,fi)
         if (obj) then
          if (fi.eq.f) then
           do j=1,3
              d=d*ten
              x(i)=y+d
              call fun(x,fi)
              if (fi.ne.f) exit
           enddo
          endif
         endif
         g(i)=(fi-f)/d
         if (obj) then
           if (center) then
              x(i)=y-d
              call fun(x,fi)
              g(i)=half*(g(i)+(f-fi)/d)
           endif
         endif
         x(i)=y
      enddo
      end

      subroutine soptions(default)

      double precision default(13)
       default(1)=-1.d0
       default(2)=1.d-4
       default(3)=1.d-6
       default(4)=15.d3
       default(5)=0.d0
       default(6)=1.d-8
       default(7)=2.5d0
       default(8)=1.d-11
       default(9)=0.d0
       default(10)=0.d0
       default(11)=0.d0
       default(12)=0.d0
       default(13)=0.d0
      end

