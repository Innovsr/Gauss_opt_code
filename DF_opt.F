c   dimension
#define   MNS   11
#define MAXVAL  10
#define MNOCC   30
#define MNBAS  700
#define   MNB  700
#define    MN  750
c   files
#define  STDIN   5
#define STDOUT   4
#define IPUNCH   7
#define   IWFN   9
#define   ENG    19
#define   RADA   17
#define   INRA   21
c   variables
#define  TINY  1.0D-10
#define  MASTER  0

c---------------------------------------------------------------------------
c         Relavistic SCF code (based on gaussian basis) 
c  This code takes three different types of basis 
c (i)    Universal Basis set (UB)
c (ii)   Well Tempered (WT)
c (iii)  Even Tempered (ET)
c  geetha 11/7/2000
c---------------------------------------------------------------------------
c             sample 'inputfile'
c  ba.dat
c  ba.out
c-------------------------------------------------------------------------
c 1st card: title
c 2nd card: name of the atom in captial
c 3rd card: total number symetry,  total s,  total p_1/2,  total p_3/2, ....
c 4th card:                         occ  s,    occ p_1/2,   occ  P_3/2,.....
c 5th card:    Kappa values
c 6th card:    J values
c 7th card:    Which method to use 
c   UB(Universal Basis)  WT (Well Tempered)  ET (Even Tempered)
c 8th card:
c  if(UB)   read alpha1,beta1
c  if(WT)   read alpha2,beta2,gamma,delta
c  if(ET)   read alpha3(isym) for each symmetry 
c           read beta3(isym)  for each symmetry
c 9th card:    charge
c 10th card:    number iteration, accuracy (give the integervalue)
c 11th card:    printing
c 12th card:   point nuclei (0) ferminuclei (1)
c 13th card:   Atomic weight and atomic number 
c-------------------------------------------------------------------------
c   Also one has to give a filename 'inputfile'
c   this contains two lines only
c 1st card:   name of the inputfile
c 2nd card:   name of the outputfile name
c-------------------------------------------------------------------------
      program main
      implicit real*8(a-h,o-z)
      character*3 state(11),chk
c     include 'mpif.h'
      common/ET1/alfa(MNS),beta(MNS),it_max
     :  /totbasis/itot,norb(200),ida(13),icom(13),idat(13)
      common/optgr2/grasp1(200),grasp2(200),grasp3(200)
     :  /che/chk
      dimension betaaa(100),bb(100),alfaaa(100),alfa_opt(15,11),
     :beta_opt(15,11),bet(10),beta_better(200),bb_best(200)
      data state/' S ',' P-',' P ',' D-',' D ',' F-',' F ',
     :                  ' G-',' G ','H-','H'/
c     call MPI_INIT(ierr)
c     START_TIME=MPI_WTIME()
c     call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
c     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      open (unit=ENG,file='grasp_E_out',status='unknown')
       idat(0)=0

       read(ENG,*)
       read(ENG,*)chk
       read(ENG,*)
       read(ENG,*)it_max,(icom(ij),ij=1,it_max)
       do i=1,it_max
       idat(i)=idat(i-1)+icom(i)
       ida(i)=idat(i)
c      print*,'idat=',i,idat(i)
       enddo

cccc total number of orbitals to be optimized, read from 'grasp_E_out' ccccccc

       iii=idat(it_max)
       itot=iii

       read(ENG,*)
       do i=1,iii
       read(ENG,*)norb(i),grasp1(i),grasp3(i),grasp2(i)
       enddo
       do i=1,itot
        if(i.le.ida(2).and.i.gt.ida(1))norb(i)=norb(i)-1
        if(i.le.ida(3).and.i.gt.ida(2))norb(i)=norb(i)-1
        if(i.le.ida(4).and.i.gt.ida(3))norb(i)=norb(i)-2
        if(i.le.ida(5).and.i.gt.ida(4))norb(i)=norb(i)-2
        if(i.le.ida(6).and.i.gt.ida(5))norb(i)=norb(i)-3
        if(i.le.ida(7).and.i.gt.ida(6))norb(i)=norb(i)-3
        if(i.le.ida(8).and.i.gt.ida(7))norb(i)=norb(i)-4
        if(i.le.ida(9).and.i.gt.ida(8))norb(i)=norb(i)-4
        if(i.le.ida(10).and.i.gt.ida(9))norb(i)=norb(i)-5
        if(i.le.ida(11).and.i.gt.ida(10))norb(i)=norb(i)-5
          enddo

cccc nn = maximum number of iteration to get ideal 'Beta value' cccc
cccc it_max = number of symmetry                                cccc
      nn=125

c     beta(1)= 1.98
c     alfa(1)= 0.0001
c     beta(2)= 1.72
c     alfa(2)= 0.0001
c     beta(3)= 1.76
c     alfa(3)= 0.0001
c     beta(4)= 1.48
c     alfa(4)= 0.0025
c     beta(5)= 1.58
c     alfa(5)= 0.0021
c     beta(6)= 1.32
c     alfa(6)= 0.0029
c     beta(7)= 1.63
c     alfa(7)= 0.0055
c     beta(8)= 1.45

cccc initial alpha and beta value for the symmetry which are not optimize cccc
      do i=1,it_max
      beta(i)= 2.73
      enddo
      do i=1,it_max
      alfa(i)= 0.0055
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc initialization    cccccccccccccccccccccccccccccccccccccccc

      b_opt_old=0
      b_opt=0
cccc  iter is the iteration generating complete sets of beta and alpha c

      iter=1

 703  b_opt_old=b_opt

cccc  it is the iteration generating beta and alpha for individual symm 

      it=1
      alfa_option=0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(alfa_option.eq.1.and.it.eq.8)goto 710
 702  print*,'beta optimisation starts for the state',state(it)
      if (it.eq.5.and.it.eq.7) goto 777

cccc  initial value for 'beta'               cccccccccccccccccccc

      beta(it)=3.00
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      print*,'iteration no',it

      do ii=2,nn
      betaaa(ii)=0.0
      enddo
 
      betaaa(1)=beta(it)
      k=1
      beta_better(1)=0

cccc   beta optimization start   cccccccccccccccccccc
      do ii=1,nn
cccc  indnt identity for beta optimization calculation cccccccccccccccc
      ident=1
ccccccccccccccccccccccccccccccccccccccccc

      print*,'beta,alfa',(beta(i),i=1,it_max),(alfa(i),i=1,it_max)
c     betaa=betaaa(ii+1)
      alfbeta=betaaa(ii)
      print*,'alfbeta',betaaa(ii)
      isk=ii
      beta_final=beta_better(k)

cccc  'test' is the main subroutine, generating Dirac Fock orbitals
      call test(b,it,ident,isk,alfbeta,indt)

cccc  'indt=' convergence criteria ccccccccccccc
      if(indt.eq.1)goto 709
ccccccccccccccccccccccccccccccccccccccccc

      bb(ii)=b
       
      if(ii.eq.1)bb_best(k)=bb(1)     
 
      if(ii.eq.1)goto 709     
     
      if(bb(ii).le.bb_best(k))then
      k=k+1    
      beta_better(k)=alfbeta
      bb_best(k)=bb(ii)
      print*,'better value of beta=',beta_better(k)
      print*,'-----------------------------------------------------'
       
      endif

cccc starting from 'beta(it=1) = 3.00' we decrease the value upto 
cccc 'beta(it=125) = 0.5' in each iteration by the 'step size = 0.02'.
cccc we store the first slope 'b'in 'bb(ii)' and update it whenever we 
cccc get a slope less than it and print the corresponding 'beta' as 
cccc 'beta_better'. The beta stored in beta_beter(k). In this way
cccc after 125 iteration we get the best 'beta' value.    

cccc  put the step size for beta cccccccccc

 709      beta(it)=beta(it)-0.02

ccccccccccccccccccccccccccccccccccccccccccc
       betaaa(ii+1)=beta(it)
        if(chk.eq.'yes')then
       goto 787
        endif
        enddo

ccccccc end of the beta optimization for the symm 'it'    ccccc
    
 708        beta(it)=beta_final
       
ccccccc alpha calculation starts cccccccccccccc
         
 710    do ii=2,nn
      bb(ii)=0.0
      alfaaa(ii)=0.0
      enddo
      bb(1)=0.0
      alfaaa(1)=alfa(it)

      print*,'alfa optimisation starts for state',state(it)

      do ii=1,nn
cccc  ident identity for alpha optimization calculation cccccccccccccccc
      ident=2
ccccccccccccccccccccccccccccccccccccccccccc
      print*,'beta,alfa',(beta(i),i=1,it_max),(alfa(i),i=1,it_max)
      alfaa=alfaaa(ii-1)
      alfbeta=alfa(ii)
      isk=ii     

cccc  'test' is the main subroutine, generating Dirac Fock orbitals
      call test(b,it,ident,isk,alfbeta,indt)

cccc  'indt=' convergence criteria ccccccccccccc
      if(indt.eq.1)goto 701

       bb(ii+1)=b
c      b_opt=bb(ii)
       if(alfaaa(ii).le.alfaaa(ii-1))goto 901

       if(ii.eq.1.or.bb(ii+1).le.bb(ii))goto 601
       if(ii.ge.3) goto 701
       if(bb(ii+1).ge.bb(ii))goto 501    
 901   if(bb(ii+1).le.bb(ii))goto 501    
       goto 701    

cccc  After getting the beta value for a specific symm. we change the alpha
cccc  value with a step size 0.0002. First we increase the value from 0.0055
cccc  and check the slope, if it is increases then start to decrease the alpha
cccc  value from 0.0055 and find the optimized alpha value at which slope 'b'is
cccc  minimum.  

 501   alfa(it)=alfa(it)-0.0002
       goto 801
 601   alfa(it)=alfa(it)+0.0002
 801   alfaaa(ii+1)=alfa(it)
       if(alfa(it).le.0) goto 701
       enddo

ccccc end of alpha optimization ccccccccccccccccccc

 701   alfa(it)=alfaa
         
       alfa_opt(iter,it)=alfa(it)
       beta_opt(iter,it)=beta(it)
       print*,'optimised value for state',state(it),'is alfa beta ='
     : ,alfa_opt(iter,it),beta_opt(iter,it)
 777   if(it.le.it_max)then
       it=it+1
       endif
       if(it.le.it_max)goto 702

cccc  a complete set of beta and alpha generation ends here cccccc

cccc  now draw a linear regression curve for complete set of optimized
cccc  beta and alpha taking all symmetry together. This is done in the 
cccc  subroutine 'all_opt' and the slope is called 'b_opt'.

       call all_opt(b_opt)
       
cccc  now slope b_opt is compared with previous slope 'b_opt_old' for
cccc  first set it is taken as '0'. If it is not the first set or we get
cccc  lesser value we continue to generate new sets till we 
cccc  get least value of slope.

       if(b_opt_old.eq.0.or.b_opt_old.gt.b_opt) iter=iter+1
       if(b_opt_old.eq.0.or.b_opt_old.gt.b_opt) goto 703
       
       print*,'the optimised alfa =',(alfa_opt(iter-1,it)
     : ,it=1,it_max)
       print*,'the corresponding beta =',(beta_opt(iter-1,it)
     :,it=1,it_max)

c      call MPI_FINALIZE(ierr)
       
 787      stop
      end
c---------------------------------------------------------
       subroutine test(b,it,ident,isk,alfbeta,indt)
c--------------------------------------------------------
      implicit real*8(a-h,o-z)
c     include 'mpif.h'

      character*2 method 
c     dimension beta3(MNS)
     
      common/quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     :      /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
     :  /ET1/alfa(MNS),beta(MNS),it_max
     :  /ET/alpha3(MNS),beta3(MNS)
     :  /totbasis/itot,norb(200),ida(13),icom(13),idat(13)
      common/optgr2/grasp1(200),grasp2(200),grasp3(200)

      data  (nkk(i),i=1,11)/1,2,2,3,3,4,4,5,5,6,6/
      data  (kap(i),i=1,11)/-1,1,-2,2,-3,3,-4,4,-5,5,-6/
      data  (unocj(i),i=1,11)/0.5d00,0.5d00,1.5d00,1.5d00,2.5d00,
     :                       2.5d00,3.5d00,3.5d00,4.5d00,4.5d00,
     :                                       5.5d00/


      call fzero
      call readinp(mocc,method)
      call setgrd
      call nucparm
      call setqic
      call nucpot

cccc 'gbasis' produce the gaussian bassis 
      call gbasis(nbasis,method)

cccc 'rinti' produce the overlap matrix and the matrix of one body part of 
cccc  the Hamiltonian which diagonalised in the first iteration in
cccc  'scfiter' and solve the equation and calculate the initial energy
cccc  and initial wave function.

      call rinti(nbasis)

cccc 'scfiter' done - 1) arrange the matrix and diagonalised it, 2)calculate
cccc  the energy and check the convergence criteria. it call subroutine 
cccc  integral and twoeint for preparing the matrixi for two body part.

      call scfiter(nbasis,indt,isk)
      if(indt.ne.1)goto 707
      if(ident.eq.1.and.indt.eq.1)goto 705
      if(ident.eq.2.and.indt.eq.1)goto 706
  705 print*,'convergence not obtained for beta',alfbeta  
      goto 707
  706 print*,'convergence not obtained for alfa',alfbeta  
  707 call output(nbasis)
      call readout(b,it)
      call daclose
      

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readout(b,it)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-h,o-z)
c     include 'mpif.h'
      dimension value(35)
      character*3 chk
c     dimension grasp1(200),grasp2(200),grasp3(200)
      common/optim/gauss1(200),gauss2(200),gauss3(200)
      common/optgr/grasp(200),gauss(200),graspp1(200),graspp2(200)
     : ,graspp3(200)
      common/optgr2/grasp1(200),grasp2(200),grasp3(200)
     :  /ET1/alfa(MNS),beta(MNS),it_max
     :  /totbasis/itot,norb(200),ida(13),icom(13),idat(13)
     :  /che/chk
 


       init=20
       jj=init+5

      if(chk.eq.'yes')then
       do i=1,itot
      print*,'grasp1',grasp1(i)
      print*,'gauss1',gauss1(i)
      print*,'grasp2',grasp2(i)
      print*,'gauss2',gauss2(i)
      print*,'grasp3',grasp3(i)
      print*,'gauss3',gauss3(i)
       enddo
      endif

c     print*,'itot',itot,ida(1),norb(1)

       
       do k=21,icom(it)+20
        gauss(k)=gauss3(k-20+idat(it-1))
       enddo

       do k=icom(it)+21,2*icom(it)+20
        gauss(k)=gauss1(k-20-icom(it)+idat(it-1))
       enddo

       do k=2*icom(it)+21,3*icom(it)+20
        gauss(k)=gauss2(k-2*icom(it)-20+idat(it-1))
       enddo

       do k=21,icom(it)+20
        grasp(k)=grasp3(k+idat(it-1)-20)
       enddo

       do k=icom(it)+21,2*icom(it)+20
        grasp(k)=grasp1(k+idat(it-1)-(20+icom(it)))
       enddo

       do k=2*icom(it)+21,3*icom(it)+20
        grasp(k)=grasp2(k+idat(it-1)-(20+2*icom(it)))
       enddo

c      do i=init+1,init+15
c      print*,'gauss',ABS(gauss(i))
c      print*,'grasp',ABS(grasp(i))
c      enddo

       do i=21,3*icom(it)+20
        value(i)=ABS(ABS(grasp(i))-ABS(gauss(i)))
        print*,'value(i)',value(i)
       enddo
       i1=0      
       do i=21,3*icom(it)+20
       i1=i1+i
       enddo
     
       value1=0.0
       
       do i=21,3*icom(it)+20
       value1=value1+value(i)
       enddo
c      print*,'value1',value1
       
       isq=0.0
       do i=21,3*icom(it)+20
       isq=isq+i*i
       enddo
c      print*,'isq',isq
       
       value2=0.0
      
       do i=21,3*icom(it)+20
       value2=value2+i*value(i)
       enddo
c      print*,'value2',value2
c      b=((i1*value1)-value2)/(i1*i1-isq)
       b=(value2)/(isq)
       print*,'b',b
       close (RADA)
       close (INRA)
       close (ENG)
       return 
       end
cc*************************************************
      subroutine all_opt(b_opt)
cc************************************************
      implicit real*8(a-h,o-z)
c     include 'mpif.h'
      dimension vall(235)
c     common/val/b
      common/optim/gauss1(200),gauss2(200),gauss3(200)
      common/optgr/graspp(200),gauss(200),graspp1(200),graspp2(200)
     : ,graspp3(200)
     :  /ET1/alfa(MNS),beta(MNS),it_max
     :  /totbasis/itot,norb(200),ida(13),icom(13),idat(13)
 
      open (unit=ENG,file='grasp_E_out',status='unknown')

       init=20
cc **************  nd=no of data points for each symmetry *************
       nd=5
cc **************  npar=no of parameter to optimised *****************
       npar=3
       iii=itot
       jj=init+nd
c      ii=(it-1)*npar*nd
c      ll=(it-1)*nd


       read(ENG,*)
       read(ENG,*)
       do i=1,iii
       read(ENG,*)graspp1(i),graspp3(i),graspp2(i)
       enddo


       do it=1,it_max
       ii=(it-1)*npar*nd
       ll=(it-1)*nd
       do k=init+1,init+iii
        gauss(k)=gauss3(k-init)
       enddo

       do k=init+1+iii,init+2*iii
        gauss(k)=gauss1(k-init-iii)
       enddo

       do k=init+1+2*iii,init+3*iii
        gauss(k)=gauss2(k-init-2*iii)
       enddo

       do k=init+1,init+iii
        graspp(k)=graspp3(k-init)
       enddo

       do k=init+1+iii,init+2*iii
        graspp(k)=graspp1(k-init-iii)
       enddo

       do k=init+1+2*iii,init+3*iii
        graspp(k)=graspp2(k-init-2*iii)
       enddo
       enddo
       do i=init+1,init+3*iii
       print*,'gauss',ABS(gauss(i))
       print*,'graspp',ABS(graspp(i))
       enddo


              

c       init=20
c       ii=it_max*15


c       do i=init+1,ii+init
c        gauss(i)=gauss(i-init)
c        grasp(i)=grasp(i-init)
c       enddo

       do i=init+1,3*iii+init
        vall(i)=ABS(ABS(graspp(i))-ABS(gauss(i)))
        print*,'vall(i)',vall(i)
       enddo
       i1=0      
       do i=init+1,3*iii+init
       i1=i1+i
       enddo
     
       vall1=0.0
       
       do i=init+1,3*iii+init
       vall1=vall1+vall(i)
       enddo
c      print*,'vall1',vall1
       
       isqr=0.0
       do i=init+1,3*iii+init
       isqr=isqr+i*i
       enddo
c      print*,'isqr',isqr
       
       vall2=0.0
      
       do i=init+1,3*iii+init
       vall2=vall2+i*vall(i)
       enddo
c      print*,'vall2',vall2
c      b_opt=((i1*vall1)-vall2)/(i1*i1-isqr)
       b_opt=(vall2)/(isqr)
       print*,'b_opt',b_opt
       close (ENG)

       return 
       end
c------------------------------------------------
      subroutine simpson(tmp,result)
c    this is the simpson integration 

      implicit real*8(a-h,o-z)
c     include 'mpif.h'
      real*8 tmp,result
      dimension tmp(MN)
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     :  /cons/zero,half,tenth,one,two,three,ten
c
      result=zero
      do l=2,n,2
      result=result+tmp(l)+tmp(l)+tmp(l+1)
      result1=result+tmp(l)+tmp(l)+tmp(l+1)
      enddo
      result=h*(result+result+tmp(1)-tmp(n))/three
      result=result+half*tmp(1)
      return
      end      

      subroutine rinti(nbasis)
c-------------------------------------------------------------------------
c     this generates overlap (sss and sll) and one-electron fock matrix,
c     where sss and sll  stands for small and large component part of the
c     overlap matrix and stored in the array emt.
c     the off-diagonal part of the fock matrix contains the kinetic energy
c     part while the diagonal part contains the nuclear contribution and
c     is stored in the array df_single. If nucopt.ne.0 then finite nucleus
c     will be considered else point nuclues.
c--------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c     include 'mpif.h'
      common/basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /overlap/sss(MNB,MNB),sll(MNB,MNB)
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /datt/nsym,nocc,nbas(MNS),ind
     :  /single/df_single(MNB,MNB)
     :  /def1/amass,z
     :  /npot/zz(MN)
     :  /skip/nskip(MNS)
     :  /ovlp/emt(MNB,MNB)
     :  /print/jprint,iout
     :  /nucop/nucopt
     :  /info2/qqf(MN,MNBAS)
     :  /info3/df_t(MNB,MNB)


      dimension atemp(MNB,MNB)
      do i=1,MNB
      do j=1,MNB
      atemp(i,j)=0.0d0
c     emt(i,j)=0.0d0
c     df_t(i,j)=0.0d0
c     df_single(i,j)=0.0
c     sss(i,j)=0.0d0
c     sll(i,j)=0.0
      enddo
      enddo
c     do i=1,n
c     zz(i)=0.0
c     enddo
c     result1=0.0
c     result2=0.0
c     result3=0.0


      spi = dsqrt(pi)
      four=two*two
      sq2=dsqrt(two)
      
      if(nucopt.eq.0)then
      do i=1,n
      zz(i)=z
      enddo
      endif

      do isym=1,nsym
      do ia=1,nbas(isym)
      ja=ia+nskip(isym)
      do ib=1,nbas(isym)
      jb=ib+nskip(isym)

      ka=ia+2*nskip(isym)
      index1=ia+nbas(isym)+2*nskip(isym)

      kb=ib+2*nskip(isym)
      index2=ib+nbas(isym)+2*nskip(isym)
c     do la = 1,n
c     tmp1(la)=0.0
c     tmp2(la)=0.0
c     tmp3(la)=0.0
c     enddo

      do la = 1,n
      tmp1(la)=gl(la,ja)*gl(la,jb)*rp(la)
      tmp2(la)=gs(la,ja)*gs(la,jb)*rp(la)
      tmp3(la)=gs(la,jb)*qqf(la,ja)*rp(la)
      enddo
      call simpson(tmp1,result1)
      call simpson(tmp2,result2)
      call simpson(tmp3,result3)
      sll(ja,jb)=result1
      sss(ja,jb)=result2
      atemp(ja,jb)=atemp(ja,jb)+c*result3
      emt(ka,kb)=result1
      emt(index1,index2) = result2
      enddo
      enddo


      do ia=1,nbas(isym)
      ja=ia+2*nskip(isym)
      ka=ia+nskip(isym)
      index1=ia+nbas(isym)+2*nskip(isym)
      do ib=1,nbas(isym)
      jb=ib+2*nskip(isym)
      kb=ib+nskip(isym)
      index2=ib+nbas(isym)+2*nskip(isym)
      df_single(index1,index2)=df_single(index1,index2)-
     :                         two*c*c*sss(ka,kb)
      enddo
      enddo

      do ia=1,nbas(isym)
      ja=ia+2*nskip(isym)
      ka=ia+nskip(isym)
      do ib=1,nbas(isym)
      kb=ib+nskip(isym)
      index2=ib+nbas(isym)+2*nskip(isym)
      df_single(ja,index2)=df_single(ja,index2)+atemp(ka,kb)
      df_single(index2,ja)=df_single(ja,index2)
      df_t(ja,index2)=df_single(ja,index2)
      df_t(index2,ja)=df_single(ja,index2)
      enddo
      enddo

c     nuclear repulsion energy

      do ia=1,nbas(isym)
      ka=ia+nskip(isym)
      ja=ia+2*nskip(isym)
      index1=ia+nbas(isym)+2*nskip(isym)
      do ib= 1,nbas(isym)
      kb=ib+nskip(isym)
      jb=ib+2*nskip(isym)
      index2=ib+nbas(isym)+2*nskip(isym)
      do la = 1,n
      tmp1(la)=gl(la,ka)*gl(la,kb)*rpor(la)*zz(la)
      tmp2(la)=gs(la,ka)*gs(la,kb)*rpor(la)*zz(la)
      enddo
      call simpson(tmp1,result1)
      call simpson(tmp2,result2)
c     matrix df_single
      df_single(ja,jb)=df_single(ja,jb)-result1
      df_single(index1,index2)=df_single(index1,index2)
     :                         -result2
      enddo
      enddo
      enddo
      mbasis=2*nbasis
      if(jprint.ne.0)then
      call matout(df_single,MNB,mbasis,mbasis,mbasis,4hH   )
      endif
      return
      end

      subroutine quad (ta,mtp,result)
c------------------------------------------------------------------------------
c   the integral in the interval zero to r(2) is computed by use of an 
c   analytical fit                                                    
c                                sigma                                
c                      f(r) = a r                                     
c   a five-point  closed  newton-cotes  formula (cf. f b hildebrand,  
c   introduction to numerical analysis, second edition, mcgraw-hill,  
c   new york, 1974, p 93)  is  used  to  compute the integral in the  
c   interval  r(2:mtp).                                              
c------------------------------------------------------------------------------
      implicit real*8(a-h, o-z)
c     include 'mpif.h'
      dimension ta(MN)
      common/cons/zero,half,tenth,one,two,three,ten
     : /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     : /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
c------------------------------------------------------------------------------
c   find first values that will permit computation of exponent
c------------------------------------------------------------------------------
c     print*,'c1,c2,c3,c4,cnc5c(2,2)',c1,c2,c3,c4,cnc5c(2,2)
      result = zero
      mtpm1 = mtp - 1
      do 3 i = 2,mtpm1
        tai = ta(i)
        if (dabs (tai) .gt. zero) then
          ip1   = i+1
          taip1 = ta(ip1)
          quott = taip1/tai
          if (quott .gt. zero) then
c------------------------------------------------------------------------------
c                          exponent from fit
c------------------------------------------------------------------------------
               frip1 = taip1/rp(ip1)
               fri   = tai  /rp(i  )
               ratio = frip1/fri
               rip1  = r (ip1)
               ri    = r (i  )
               sigma = dlog (ratio)/dlog (rip1/ri)
c------------------------------------------------------------------------------
c       analytical integration and error estimate for interval r(1:i)
c------------------------------------------------------------------------------
               fri    = ri*fri
               result = fri/(sigma+one)
c------------------------------------------------------------------------------
c                       set the tail to zero
c------------------------------------------------------------------------------
               do 1 loc = 1,3
                  ta(mtp+loc) = zero
    1          continue
c------------------------------------------------------------------------------
c              newton-cotes quadature for the remainder
c------------------------------------------------------------------------------
               result = result+c1*tai
               do 2 loc = ip1,mtp,4
                  result = result+c2*(ta(loc  )+ta(loc+2))
     :                           +c3* ta(loc+1)
     :                           +c4* ta(loc+3)
    2          continue
               if (mod (mtp-i,4) .eq. 0) result = result-c1*ta(mtp)
               goto 4
            endif
         endif
    3 continue
c------------------------------------------------------------------------------
c          no value which will permit computation of exponent
c------------------------------------------------------------------------------
      result = zero
      do j = 2,5
         do i = 2,5
            c5num(i,j) = 0.0
        enddo
      enddo
    4 return
      end
      
      subroutine factt(gam)
c---------------------------------------------------------------------
c calculates the logs  of factorials required by the racah coefficient 
c routine  dracah. written by n.s. scott.
c---------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c     include 'mpif.h'
      dimension gam(500)
      common /cons/zero,half,tenth,one,two,three,ten
      data thirty/3.0d 01/
*
      do i=1,500
       gam(i)=0.0
      enddo

      gam(1) = one
      gam(2) = one
      x = two
      do i = 3,30
         gam(i) = gam(i-1)*x
         x = x+one
      enddo
      do  i = 1,30
         gam(i) = dlog(gam(i))
      enddo
      x = thirty
      do i = 31,500
         gam(i) = gam(i-1)+dlog(x)
         x = x+one
      enddo
      return
*
      end
      subroutine init(ii) 
c----------------------------------------------------------------------
c   this subroutine will initialiaze the arrays and the constants in  *
c   dummy blocks dumm and datt                                        *
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c     include 'mpif.h'
      common /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /datt/nsym,nocc,nbas(MNS),ind
      do 10 ia=1,MN
      tmp1(ia)=zero
      tmp2(ia)=zero
      tmp3(ia)=zero
      tmp4(ia)=zero
 10   continue
      cn2=zero
      cn3=zero
      cn4=zero
      return
      end
      subroutine scf(nbasis)
c------------------------------------------------------------------
c     This performs the scf iteration.
c------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c     include 'mpif.h'
      logical set
      character*3 nh_i
      real*8 df_mat,swap
      dimension can(MNB,MNB)
      dimension  drs(MNB,MNB),sa(MNB),df_mat(MNB,MNB)
      common/basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /overlap/sss(MNB,MNB),sll(MNB,MNB)
     :  /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
     :  /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     :  /diag/eigv(MNB,MNB),eng(MNB)
     :  /datt/nsym,nocc,nbas(MNS),ind
     :  /single/df_single(MNB,MNB)
     :  /double/df_two(MNB,MNB)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /conv/crit,maxit,npower
     :  /print/jprint,iout
     :  /jv/occ(MNOCC),orbe(MNB),iq(MNB)
     :  /skip/nskip(MNS)
     :  /dens/p(MNB,MNB)
     :  /ovlp/emt(MNB,MNB)
     :  /info1/occk(MNOCC,MNS)
     :  /info/nocorb(MNOCC),nucorb(MNB)
     :  /vector/vec(MNB,MNB),oe(MNB,MNS)
     :  /info3/df_t(MNB,MNB)

      dimension a(MNB,MNB),st(MNB,MNB),b(MNB,MNB)
      dimension pold(MNB,MNB),stemp(MNB,MNB)
      dimension nh_i(MNS),scr(MNB,MNB)

      data nh_i/' S ',' P-',' P ',' D-',' D ',' F-',' F ',
     :                  ' G-',' G ','H-','H'/

c--------------------------------------------------------------------
c       initialize arrays and constants required for the computation.
c--------------------------------------------------------------------
      tcsqr=-two*c*c
      set = .true.
      iter=0
      en=zero
      eold=zero
      delta=zero


      do ia=1,MNB
      sa(ia)=zero
      do ib=1,MNB
      can(ib,ia)=zero
      drs(ib,ia)=zero
      p(ia,ib)=zero
      enddo
      enddo

      mbasis=2*nbasis

      do isym=1,nsym
      ndim=nbas(isym)
      ndim2=2*ndim
 
      do ia=1,ndim2
      ja=ia+2*nskip(isym)
      do ib=1,ndim2
      jb=ib+2*nskip(isym)
      stemp(ia,ib)=emt(ja,jb)
      enddo
      enddo
     
      if(jprint.ne.0)call matout(stemp,MNB,ndim2,ndim2,ndim2,4hS   )

      do i=1,ndim2
      do j=1,ndim2
      drs(i,j)=stemp(i,j)
      enddo
      enddo

      call eigen (drs,scr,ndim2,sa,ndim2,bnorm,MNB)
c-------------------------------------------------------------------
c       canonical transformation
c------------------------------------------------------------------
      do ia=1,ndim2
      do ib=1,ndim2
      if(dabs(sa(ib)).gt.TINY) then
      can(ia,ib)=drs(ia,ib)/dsqrt(dabs(sa(ib)))
      st(ib,ia)=can(ia,ib)
      else
      endif
      enddo
      enddo

      if(jprint.ne.0)call matout(st,MNB,ndim2,ndim2,ndim2,4hXT  )

      if(jprint.ne.0)
     :          call matout(df_single,MNB,mbasis,mbasis,mbasis,4hH   )

      call twoe(isym,nbasis,iter)

      do ic=1,ndim2
      jc=ic+2*nskip(isym)
      do id=1,ndim2
      jd=id+2*nskip(isym)
      a(ic,id)=df_single(jc,jd)+df_two(jc,jd)
      enddo
      enddo

      if(jprint.ne.0)call matout(a,MNB,ndim2,ndim2,ndim2,4hFOCK )

      if(jprint.ne.0)
     :call matout(df_two,MNB,mbasis,mbasis,mbasis,4hVHF )

      call mult(a,can,b,MNB,ndim2)
      call mult(st,b,df_mat,MNB,ndim2)
      if(jprint.ne.0)
     :      call matout(df_mat,MNB,ndim2,ndim2,ndim2,6hFRIME )
      call eigen (df_mat,scr,ndim2,eng,ndim2,bnorm,MNB)
      do i=1,ndim2
      do j=1,ndim2
      drs(i,j)=df_mat(i,j)
      enddo
      enddo


      if(jprint.ne.0)
     : call matout(drs,MNB,ndim2,ndim2,ndim2,4hC   )
      call mult(can,drs,eigv,MNB,ndim2)
  
      loc = 1
      last=ndim-1
      dowhile(loc.le.last)
      samin = eng(loc)
      locus = loc
      do  ia = loc,ndim
      if (eng(ia).lt.samin) then
      samin = eng(ia)
      locus = ia
      endif
      enddo
      if (locus .gt. loc) then
      swap = eng(loc)
      eng(loc) = eng(locus)
      eng(locus) = swap
      do  ib = 1,ndim2
      swap = eigv(ib,loc)
      eigv(ib,loc) = eigv(ib,locus)
      eigv(ib,locus) = swap
      enddo
      endif
      loc = loc + 1
      enddo
c-----------------------------------------------------------------
      do ia=1,ndim2
      do ib=1,ndim2
      b(ia,ib)=zero
      enddo
      b(ia,ia)=eng(ia)
      oe(ia,isym)=eng(ia)
      enddo

      if(jprint.ne.0)call matout(b,MNB,ndim2,ndim2,ndim2,4hE   )

      if(jprint.ne.0)call matout(eigv,MNB,ndim2,ndim2,ndim2,4hCP  )

      do ia=1,ndim2
      ja=ia+2*nskip(isym)
      do ib=1,ndim2
      jb=ib+2*nskip(isym)
      pold(ja,jb)=p(ja,jb)
      p(ja,jb)=0.0d0
      nocc=nocorb(isym)
      if(nocc.ne.0)then
      do ibas=1,nocc
      twojp1=two*occ(isym)+one
      p(ja,jb)=p(ja,jb)+twojp1*eigv(ia,ibas)*eigv(ib,ibas)
      enddo
      endif
      vec(ja,jb)=eigv(ia,ib)
      enddo
      enddo

      enddo

      if(jprint.ne.0)call matout(p,MNB,mbasis,mbasis,mbasis,4hP   )
      
      call integral(nbasis)
      eold=en
      en=0.0d0
      do i=1,mbasis
      do j=1,mbasis
      en=en+0.5d0*p(i,j)*(two*df_single(i,j)+df_two(i,j))
      print*,'**********'
      enddo
      enddo
      
      write(STDOUT,204)en    
      do jsym=1,nsym
      do ibas=1,nbas(jsym)
      write(STDOUT,206)ibas,nh_i(jsym),oe(ibas,jsym),
     :                 oe(ibas+nbas(jsym),jsym)
      

  
      enddo
      enddo

      return
 199  format(/,1x,'CONVERGENCE DATA',/,
     1        ' MAXIMUM NO. OF ITERATION=',i6,/,
     2        ' CONVERGENCE CRITERION   = 1.0D-',i2,//)
 200  format(/,1x,'CYCLE      DENSITY CONV',/)
 201  format(i3,1x,(2x,d20.7))
 202  format(///,4x,33hSCF FAILS TO CONVERGES AT CYCLE  , i4)
 203  format(///,4x,22hSCF CONVERGES AT CYCLE, i4)
 204  format(///,4x,20hELECTRONIC ENERGY=  ,d20.11)
 205  format(///,4x,32h ORBITAL ENERGIES (+ve and -ve) ,i4,///)
 206  format(/,4x,I2,A,2D20.11)
      end            
c---------------------------------------------------------------------
      subroutine output(nbasis)
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c     include 'mpif.h'
      common/wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     : /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     : /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     : /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     : /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     : /datt/nsym,nocc,nbas(MNS),ind
     : /vector/eigv(MNB,MNB),eng(MNB,MNS)
     :  /jv/occ(MNOCC),orbe(MNB),iq(MNB)
     :  /skip/nskip(MNS)
     :  /print/jprint,iout
     :  /kapp/kappa(MNB),kappc(MNB),iqc(MNOCC)
     : /freez/enmin(MNS),enmax(MNS)
      common/optim/gauss1(200),gauss2(200),gauss3(200)
     :  /totbasis/itot,norb(200),ida(13),icom(13),idat(13)

      dimension ns(11)
 
      write(STDOUT,208)
c-----------------------------------------------------------------------
c            now construct the gaussian orbitals on a grid.
c-----------------------------------------------------------------------
c      print*,'pf,qf',pf(3,100),qf(3,100)
      write(IWFN)h,n
      write(IWFN)(r(i),i = 1,n),(rp(i),i = 1,n),(rpor(i),i=1,n)
      do isym=1,nsym

c     write only the relevant orbitals i.e skip writting the
c     negavtive energy orbitals
******************************************
*   if you like to have expansion coeff
*****************************************
c         write(IWFN)nbas(isym),nskip(isym)
*****************************************
      do  ia=1,nbas(isym)
c   initialise 
               do ic=1,n
                       pf(ic,ia)=0.0d0
                       qf(ic,ia)=0.0d0
               enddo
       enddo
       do  ia=1,nbas(isym)
       ka=ia+nskip(isym)
       ja=ia+2*nskip(isym)
 	do ib =1,nbas(isym)
        jb=ib+2*nskip(isym)
         kb=ib+nskip(isym)
          do  ic=1,n
            
            pf(ic,ia)=pf(ic,ia)+eigv(jb,ja)*gl(ic,kb)
            qf(ic,ia)=qf(ic,ia)+eigv(jb+nbas(isym),ja)*gs(ic,kb)
	  enddo
******************************************
*   if you like to have expansion coeff
c         write(IWFN)eigv(jb,ja),eigv(jb+nbas(isym),ja)
*******************************************          
        enddo


c     if (iout.eq.0) then
      write(15,*)ia,eng(ia,isym)
      write(IWFN)eng(ia,isym)
      write(IWFN)(pf(ii,ia),ii=1,n),(qf(ii,ia),ii=1,n)
      do ib =1,ia
      kb=ib+nskip(isym)
      result = zero
      do ic = 1,n
      tmp1(ic) = (pf(ic,ib)*pf(ic,ia)+qf(ic,ib)*
     :                     qf(ic,ia))*rp(ic)
      enddo
c  integrate tmp1 for orthogonality
      call quad(tmp1,n,result)
      r0=rint(ia,ib,0)
      r1=0
      r2=0
      r3=0
      if(ia.eq.ib)r1=rint(ia,ib,1)
      if(ia.eq.ib)r2=rint(ia,ib,2)
      if(ia.eq.ib)r3=rint(ia,ib,-1)
      write(STDOUT,'(2i4,2x,4(2x,D15.9))')ka,kb,r0,r1,r2,r3
      enddo

c     write(15,*)ia,eng(ia,isym),enmax(isym)

      nbas2=0
      nbas3=0
      nbas4=0
      nbas5=0
      nbas6=0
      nbas7=0
      nbas8=0
      nbas9=0
      nbas10=0
      do i=1,2
      nbas2=nbas2+nbas(i)
      enddo
      do i=1,3
      nbas3=nbas3+nbas(i)
      enddo
      do i=1,4
      nbas4=nbas4+nbas(i)
      enddo
      do i=1,5
      nbas5=nbas5+nbas(i)
      enddo
      do i=1,6
      nbas6=nbas6+nbas(i)
      enddo
      do i=1,7
      nbas7=nbas7+nbas(i)
      enddo
      do i=1,8
      nbas8=nbas8+nbas(i)
      enddo
      do i=1,9
      nbas9=nbas9+nbas(i)
      enddo
      do i=1,10
      nbas10=nbas10+nbas(i)
      enddo
      ns(1)=0
      ns(2)=nbas(1)
      ns(3)=nbas2
      ns(4)=nbas3
      ns(5)=nbas4
      ns(6)=nbas5
      ns(7)=nbas6
      ns(8)=nbas7
      ns(9)=nbas8
      ns(10)=nbas9
      ns(11)=nbas10
cc          <r>  reading for optimisation
ccc    S    symmetry
      do i=1,itot
        if(i.le.ida(1))nsl=1
        if(i.le.ida(2).and.i.gt.ida(1))nsl=2
        if(i.le.ida(3).and.i.gt.ida(2))nsl=3
        if(i.le.ida(4).and.i.gt.ida(3))nsl=4
        if(i.le.ida(5).and.i.gt.ida(4))nsl=5
        if(i.le.ida(6).and.i.gt.ida(5))nsl=6
        if(i.le.ida(7).and.i.gt.ida(6))nsl=7
        if(i.le.ida(8).and.i.gt.ida(7))nsl=8
        if(i.le.ida(9).and.i.gt.ida(8))nsl=9
        if(i.le.ida(10).and.i.gt.ida(9))nsl=10
        if(i.le.ida(11).and.i.gt.ida(10))nsl=11
      if(ka.eq.norb(i)+ns(nsl).and.kb.eq.norb(i)+ns(nsl)) gauss2(i)=r1
      enddo

cc          <1/r>  reading for optimisation

ccc    S    symmetry
      do i=1,itot
        if(i.le.ida(1))nsl=1
        if(i.le.ida(2).and.i.gt.ida(1))nsl=2
        if(i.le.ida(3).and.i.gt.ida(2))nsl=3
        if(i.le.ida(4).and.i.gt.ida(3))nsl=4
        if(i.le.ida(5).and.i.gt.ida(4))nsl=5
        if(i.le.ida(6).and.i.gt.ida(5))nsl=6
        if(i.le.ida(7).and.i.gt.ida(6))nsl=7
        if(i.le.ida(8).and.i.gt.ida(7))nsl=8
        if(i.le.ida(9).and.i.gt.ida(8))nsl=9
        if(i.le.ida(10).and.i.gt.ida(9))nsl=10
        if(i.le.ida(11).and.i.gt.ida(10))nsl=11
      if(ka.eq.norb(i)+ns(nsl).and.kb.eq.norb(i)+ns(nsl)) gauss3(i)=r3
      enddo
c     enddo
c     endif
 8888  continue
c     endif
      enddo
      enddo

c     rewind(IWFN)
c     read(IWFN)h,n
c     read(IWFN)(r(i),i = 1,n),(rp(i),i = 1,n),(rpor(i),i=1,n)
c     do i=1,nbasis
c     read(IWFN)orben
c     read(IWFN)(pf(ii,i),ii=1,n),(qf(ii,i),ii=1,n)
c     enddo

      return
 208  format(2x,'THESE ARE THE ORTHOGONALITY OF THE  ',
     :  1x,'ORBITALS',//,2x,'Sl. No.','    <A|B>      ',
     :  '    <A|R|B>     ','    <A|R*R|B>   ',
     :  '     <A|(1/R)|B>',//)
 210  format(///,2x,'THESE ARE THE DIPOLE MOMENTS OF THE  ',
     :  1x,'ORBITALS',//,2x,'Sl. No.',5x,'DIPOLE MOMENTS'/)
      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real*8 function clrx (kappaa,k,kappab)
c----------------------------------------------------------------------
c   the value of clrx is the 3-j symbol:
c
c                    ( ja        k        jb  )
c                    ( 1/2       0       -1/2 )
c 
c the  k's are kappa angular quantum numbers. the formula is taken from d m
c brink and g r satchler, <angular momentum>, second edition (oxford: 
c clarendon press, 1968), p 138.   the logarithms of the first  mfact  
c factorials must be available in  common/facts/ for this program to function 
c correctly. note that  n!  is stored in fact(n+1).
c----------------------------------------------------------------------
      implicit real*8(a-h, o-z)
c     include 'mpif.h'
      dimension gam(500)
      common/cons/zero,half,tenth,one,two,three,ten

       call factt(gam) 
c----------------------------------------------------------------------
c               determine the absolute values of the kappas
c----------------------------------------------------------------------
      ka =iabs (kappaa)
      kb =iabs (kappab)
c----------------------------------------------------------------------
c                   perform the triangularity check
c----------------------------------------------------------------------
      if ((iabs(ka-kb) .le. k) .and. (ka+kb-1 .ge. k)) then
c----------------------------------------------------------------------
c   triangularity satisfied; compute the 3j coefficient
c   begin with the logarithm of the square of the leading term
c----------------------------------------------------------------------
         exptrm = -dlog (dble (ka*kb))
c----------------------------------------------------------------------
c    compute the logarithm of the square root of the leading term and the 
c    factorial part that doesn't depend on the parity of ka+kb+k (the delta 
c    factor)
c----------------------------------------------------------------------
         kapkb = ka+kb
         kabkp = kapkb+k
         kamkb = ka-kb
         kbmka = kb-ka
         exptrm = half*(exptrm+gam(kapkb-k  )+gam(kamkb+k+1)
     :                        +gam(kbmka+k+1)-gam(kabkp  +1) )
c----------------------------------------------------------------------
c   the remainder depends on the parity of ka+kb+k
c----------------------------------------------------------------------
         if (mod (kabkp,2) .eq. 0) then
c----------------------------------------------------------------------
c   computation for even parity case
c   include the phase factor: a minus sign if necessary
c----------------------------------------------------------------------
            if (mod (3*kabkp/2,2) .eq. 0) then
               clrx =  one
            else
               clrx = -one
            endif
c----------------------------------------------------------------------
c   include the contribution from the factorials
c----------------------------------------------------------------------
            exptrm = exptrm+gam((kabkp  +2)/2)-gam((kapkb-k  )/2)
     :                     -gam((kamkb+k+2)/2)-gam((kbmka+k+2)/2)
*
         else
c----------------------------------------------------------------------
c                   computation for odd parity case
c         include the phase factor: a minus sign if necessary
c----------------------------------------------------------------------
            if (mod ((3*kabkp-1)/2,2) .eq. 0) then
               clrx =  one
            else
               clrx = -one
            endif
c----------------------------------------------------------------------
c   include the contribution from the factorials
c----------------------------------------------------------------------
            exptrm = exptrm+gam((kabkp  +1)/2)-gam((kapkb-k+1)/2)
     :                     -gam((kamkb+k+1)/2)-gam((kbmka+k+1)/2)
*
         endif
c----------------------------------------------------------------------
c   final assembly
c----------------------------------------------------------------------
         clrx = clrx*dexp (exptrm)
*
      else
c----------------------------------------------------------------------
c   triangularity violated; set the coefficient to zero
c----------------------------------------------------------------------
         clrx = zero
*
      endif
*
      return
      end
      subroutine setqic
c--------------------------------------------------------------------------
c this  subroutine sets up the coefficients for subroutines 
c quad, rinti, yzk 
c--------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c     include 'mpif.h'
      dimension c5num1(5,5)
      common/cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
      common/quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)

c------------------------------------------------------------------------
c   five-point newton-cotes coefficients for closed integration. 
c   expressed as  rational numbers
c-------------------------------------------------------------------------
      data (c5num1(i,2),i = 1,5)/251.0d 00,646.0d 00,
     :       -264.0d 00,106.0d 00,-19.0d 00/
      data (c5num1(i,3),i = 1,5)/232.0d 00,992.0d 00,
     :        192.0d 00,32.0d 00,-8.0d 00/
      data (c5num1(i,4),i = 1,5)/243.0d 00,918.0d 00,
     :        648.0d 00,378.0d 00,-27.0d 00/
      data (c5num1(i,5),i = 1,5)/224.0d 00,1024.0d 00,
     :        384.0d 00,1024.0d 00,224.0d 00/
      data c5den/ 720.0d 00/


c------------------------------------------------------------------------
c                  newton-cotes coefficients for yzk and quad
c-----------------------------------------------------------------------
c      c1=0.0
c      c2=0.0
c      c3=0.0
c      c4=0.0
c      factor=0.0
      factor = h/c5den

      do 5 j = 2,4
         do 4 i = 2,5
c          print*,'c5num(i,j)',i,j,c5num1(i,j)
            cnc5c(i,j) = factor*c5num1(i,j)
    4    continue
    5 continue
*
      c1 = factor*c5num1(1,5)
      c2 = factor*c5num1(2,5)
      c3 = factor*c5num1(3,5)
      c4 = c1+c1
*
      do 7 j = 2,5
         do 6 i = 2,5
            c5num(i,j) = c5num1(i,j)/c5den
c          print*,'c5num(i,j) new ',i,j,c5num(i,j),c5num(i,j)
    6    continue
    7 continue
      pi = 3.141592653589790d 000
      c = 1.3703598950000d 002
c     c = 1.00000d006
      zero=0.0d 00
      half=0.5d 00 
      tenth=0.1d 00
      one=1.0d 00
      two=2.0d 00
      three=3.0d 00
      ten=10.0d 00
      return
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer function idbf(n,i)
c     calculates double factorial
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c     include 'mpif.h'
      if(n.lt.-1) then
      write(*,*)'n.lt.-1',n
      stop
      else
        if((n.eq.0).or.(n.eq.-1)) then
          idbf=1
        else
          idbf=1
          do j=n,1,-i
            idbf=j*idbf
          enddo
        endif
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine nucpot
c----------------------------------------------------------------------
c evaluate the nuclear potential for point and fermi models.
c subroutines called: es.
c-----------------------------------------------------------------------
      implicit doubleprecision (a-h, o-z)
c     include 'mpif.h'
      logical set
*
      common /def1/amass,z
     :  /npar/parm(2),nparm
     :  /npot/zz(MN)
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
      common/quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)


      data pi/3.141592653589793d0/

c-----------------------------------------------------------------------
c                                 point nucleus
c-----------------------------------------------------------------------
      if (nparm.eq.0) then
         do 1 i = 1,MN
            zz(i) = z
    1    continue
c-----------------------------------------------------------------------
c                               fermi distribution
c-----------------------------------------------------------------------
      elseif (nparm.eq.2) then
         c = parm(1)
         a = parm(2)
         abc = a/c
         tabc = two*abc
         abc2 = abc*abc
         thabc2 = three*abc2
         abc3 = abc2*abc
         cba = c/a
         pi2 = pi*pi
         hpiac2 = half*pi2*abc2
         six = two*three
         h3 = half*three
         h3php = h3+hpiac2
         call es (-cba,s2mcba,s3mcba)
         sabc3 = six*abc3
         dmsas = -sabc3*s3mcba
         en = one + abc2*pi2 + dmsas
         zbn = z/en
*
         do 2 i = 1,MN
            ri = r(i)
            rmc = ri-c
            rmcba = rmc/a
            rbc = ri/c
            if (rbc.le.one) then
               call es (rmcba,s2rcba,s3rcba)
               zz(i) = zbn*( dmsas + sabc3*s3rcba
     :                      +rbc*( h3php-thabc2*s2rcba
     :                            -half*rbc*rbc) )
            else
               if (.not. set) then
                  nnuc = i
                  set = .true.
             endif
               call es (-rmcba,s2rcba,s3rcba)
               zz(i) = z * ( one
     :                        +thabc2 * ( rbc *s2rcba
     :                                   +tabc*s3rcba ) / en )
            endif
    2    continue
      endif
      return
      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine es (f,s2f,s3f)
c------------------------------------------------------------------------
c                        evaluate the sum of the series
c
c          k(f)=sum_(0)^(infinity)(-1)^n  exp(n*f)/ (n^k)
c
c       for k = 2, 3 to machine precision. this is a utility 
c------------------------------------------------------------------------
      implicit doubleprecision (a-h, o-z)
c     include 'mpif.h'
      common /cons/zero,half,tenth,one,two,three,ten
*
      n = 0
      s2f = zero
      s3f = zero
      fase = one
    1 n = n+1
         en = dble (n)
         obn = one/en
         fase = -fase
         enf = dexp (en*f)
         term2 = fase*enf*obn*obn
         term3 = term2*obn
         s2last = s2f
         s2f = s2f+term2
         s3f = s3f+term3
         if (dabs (s2f) .ne. dabs (s2last)) goto 1
      return
*
      end
c------------------------------------------------------------
      subroutine mult(a,b,c,im,m)
c    matrice multiplication
c------------------------------------------------------------
      implicit double precision (a-h, o-z)
c     include 'mpif.h'
      dimension a(im,im),b(im,im),c(im,im)
    
      do i=1,im
        do j=1,im
        c(i,j)=zero
        enddo
       enddo
      do i=1,m
          do j=1,m
c          c(i,j)=0.0d0
            do k=1,m
            c(i,j)=c(i,j)+a(i,k)*b(k,j)
            end do
          end do
       end do
       return
       end
c------------------------------------------------------------
      subroutine fzero
c    initialise
c------------------------------------------------------------
      implicit real*8(a-h,o-z)
c     include 'mpif.h'
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
     :  /overlap/sss(MNB,MNB),sll(MNB,MNB)
     :  /datt/nsym,nocc,nbas(MNS),ind
     :  /single/df_single(MNB,MNB)
     :  /npar/parm(2),nparm
      common /def1/amass,z
     :  /npot/zz(MN)
     :  /jv/occ(MNOCC),orbe(MNB),iq(MNB)
     :  /totbasis/itot,norb(200),ida(13),icom(13),idat(13)

      do i=1,MN
                 r(i)=zero
                 rp(i)=zero
                 rpor(i)=zero
                 tmp1(i)=zero
                 tmp2(i)=zero
                 tmp3(i)=zero
                 tmp4(i)=zero
                 zz(i)=zero
           do j=1,MNBAS
                 pf(i,j)=zero
                 qf(i,j)=zero
                 gl(i,j)=zero
                 gs(i,j)=zero
           enddo
      enddo
      do i=1,MNS
                 nsy(i)=0
                 unocl(i)=zero
                 unocj(i)=zero
                 occ(i)=zero
         do j=1,MNBAS
                 alpha(j,i)=zero
         enddo
      enddo
      do i=1,MNB
         do j=1,MNB
                 sss(i,j)=zero
                 sll(i,j)=zero
                 df_single(i,j)=zero
         enddo
      enddo
      return
      end
c***********************************************************************
      subroutine slater(ia,ib,ic,id,k,rkll,rkss,rkls,rksl)

*   the value of this  function is the  slater integral  as normally   *
*   defined in terms of the four sets of quantum numbers  a,b,c,d.     *
*                                                                      *
*   subroutines called: quad, yzk.                                     *
*                                                                      *
c***********************************************************************
*
      implicit real*8(a-h, o-z)
c     include 'mpif.h'
c     dimension rhop(MN),rttk(MN),wk(MN),temp(MN),zk(MN)
      common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)

      mtp=n
      do i=1,mtp
         tmp1(i)=0
         tmp2(i)=0
         tmp3(i)=0
         tmp4(i)=0
      enddo

c     call yzk(k,ib,id)
      do i = 2,mtp
         tmp3(i)=gl(i,ia)*gl(i,ic)*rpor(i)*tmp1(i)
      enddo
      call quad (tmp3,mtp,rkll)

      do i = 2,mtp
         tmp3(i)=gs(i,ia)*gs(i,ic)*rpor(i)*tmp2(i)
      enddo

      call quad (tmp3,mtp,rkss)

      do i = 2,mtp
         tmp3(i)=gs(i,ia)*gs(i,ic)*rpor(i)*tmp1(i)
      enddo

      call quad (tmp3,mtp,rksl)

      do i = 2,mtp
         tmp4(i)=gl(i,ia)*gl(i,ic)*rpor(i)*tmp2(i)
      enddo

      call quad (tmp4,mtp,rkls)

      return
      end

      subroutine matout(a,im,in,m,n,label)
      implicit double precision (a-h,o-z)
c     include 'mpif.h'
      character*6 label
      dimension a(im,im)
 
c       do i=1,MNB
c         do j=1,MNB
c       a(i,j)=zero
c         enddo
c       enddo
      ihigh=0
      low=0
  10  low=ihigh+1
      ihigh=ihigh+5
      ihigh=min0(ihigh,n)
              write(STDOUT,20)label,(i,i=low,ihigh)
  20  format(///,3x,5h THE ,a6,6h ARRAY,/15x,5(10x,i3,6x)//)
      do i=1,m
              write(STDOUT,40) i,(a(i,j),j=low,ihigh)
      enddo
  40  format(i8,5x,8(1x,d12.6))
      if(n-ihigh)50,50,10
  50  return
      end

      subroutine scfout
      write(STDOUT,9990)
      write(STDOUT,9999)
      write(STDOUT,9998)
      write(STDOUT,9998)
      write(STDOUT,9998)
      write(STDOUT,9997)
      write(STDOUT,9998)
      write(STDOUT,9996)
      write(STDOUT,9998)
      write(STDOUT,9995)
      write(STDOUT,9998)
      write(STDOUT,9998)
      write(STDOUT,9998)
      write(STDOUT,9999)

9999  format(20x,'**************************************')
9998  format(20x,'*                                    *')
9997  format(20x,'*        Relativistic SCF Run        *')
9996  format(20x,'*            written by              *')
9995  format(20x,'*   P. K.  Panda and R. K. Chaudhuri *')
9990  format(//)
      return
      end
c*****************************************************************
c     generate the basis functions on a grid point. It generates
c     (at present) geometrical basis but can be easily transformed
c     to a general one and also to contracted gaussian type.
c*****************************************************************
      subroutine gbasis(nbasis,method)
      implicit real*8(a-h,o-z)
c     include 'mpif.h'
      real*8 norm_large,norm_small
      character*2 method
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
     :  /overlap/sss(MNB,MNB),sll(MNB,MNB)
     :  /datt/nsym,nocc,nbas(MNS),ind
     :  /single/df_single(MNB,MNB)
     :  /npar/parm(2),nparm
      common /def1/amass,z
     :  /npot/zz(MN)
     :  /jvalue/orbj(MNOCC,MNS),iocc(MNBAS)
     :  /skip/nskip(MNS)
     :  /kapp/kappa(MNB),kappc(MNB),iqc(MNOCC)
     :  /info1/occk(MNOCC,MNS)
     :  /jv/occ(MNOCC),orbe(MNB),iq(MNB)
     :  /info/nocorb(MNOCC),nucorb(MNB)
     :  /info2/qqf(MN,MNBAS)
c    :  /UB/alpha1,beta1
c    :  /WT/alpha2,beta2,gamma,delta
     :  /ET/alpha3(MNS),beta3(MNS)



      dimension rnk(MN)
      spi=dsqrt(pi)
      four = two*two
      onehalf = one + half
      nbasis=0
      nskip(1)=0
      do isym=1,nsym
      nbasis=nbasis+nbas(isym)
      nskip(isym+1)=nskip(isym)+nbas(isym)
      enddo
  
      ii=0
      do isym=1,nsym
      do ibas=1,nbas(isym)
      ii=ii+1
      kappa(ii)=nak(isym)
      orbe(ii)=occ(isym)
      if(kappa(ii).eq.-1)iq(ii)=1
      if(kappa(ii).eq. 1)iq(ii)=-1
      if(kappa(ii).eq.-2)iq(ii)= 1
      if(kappa(ii).eq. 2)iq(ii)=-1
      if(kappa(ii).eq.-3)iq(ii)= 1
      if(kappa(ii).eq. 3)iq(ii)=-1
      if(kappa(ii).eq.-4)iq(ii)= 1
      if(kappa(ii).eq. 4)iq(ii)=-1
      if(kappa(ii).eq.-5)iq(ii)= 1
      if(kappa(ii).eq. 5)iq(ii)=-1
      if(kappa(ii).eq.-6)iq(ii)= 1
      enddo
      enddo
 
      ii=0
      do isym=1,nsym
      nn=nocorb(isym)
      if(nn.ne.0)then
      do ibas=1,nocorb(isym)
      ii=ii+1
      occl(ii)=occ(isym)
      occk(ibas,isym)=occ(isym)
      kappc(ii)=nak(isym)
      if(kappc(ii).eq.-1)iqc(ii)=1
      if(kappc(ii).eq. 1)iqc(ii)=-1
      if(kappc(ii).eq.-2)iqc(ii)= 1
      if(kappc(ii).eq. 2)iqc(ii)=-1
      if(kappc(ii).eq.-3)iqc(ii)= 1
      if(kappc(ii).eq. 3)iqc(ii)=-1
      if(kappc(ii).eq.-4)iqc(ii)= 1
      if(kappc(ii).eq. 4)iqc(ii)=-1
      if(kappc(ii).eq.-5)iqc(ii)= 1
      if(kappc(ii).eq. 5)iqc(ii)=-1
      if(kappc(ii).eq.-6)iqc(ii)= 1
      enddo
      endif
      enddo

      nocc=ii
             write(STDOUT,200)nbasis
             write(STDOUT,201)
      do isym=1,nsym
             write(STDOUT,202)isym,nbas(isym),nskip(isym)
      enddo

      do i=1,nbasis
             mf(i)=n
             mg(i)=n
      enddo

      if(method.eq.'UB')then
      iopt=1
      endif
      if(method.eq.'WT')then
      iopt=2
      endif
      if(method.eq.'ET')then
      iopt=3
      endif
 
      sq2=dsqrt(two)
      call init(1)
c-----------------------------------------------------------------
c              generate the gaussians on a grid
c         one thing to be noted: these are normalized  guassians. 
c-----------------------------------------------------------------
c  here we give different alpha3 for different symmetries. So alpha3
c  is also given thro  the common block instead thro arguments.
      do isym=1,nsym

c----------------------------------------------------------------
c  here the value of the exponential factor for the 
c                 gaussians are computed.
c----------------------------------------------------------------

      do ib=1,n
      rnk(ib) = r(ib)**nkk(isym)
      enddo
      do ic=1,nbas(isym)
      jc=ic+nskip(isym)

C here the alpha is calculated using alpha(i) and beta(i) for the 
C three different types of basis given by UB,WT and ET.

      goto (1,2,3),iopt
   1  alpha(jc,isym) = alpha1*(beta1**(ic-1))
      goto 4

   2  alpha(jc,isym) = alpha2*(beta2**(ic-1))*(1+ gamma *
     :               (ic/nbas(isym))**delta)
      goto 4

   3  alpha(jc,isym) = alpha3(isym)*(beta3(isym)**(ic-1))
   4  continue
      enddo
      dkap=dble(kap(isym))
      dnkk=dble(nkk(isym))
      ifact=idbf(2*nkk(isym)-1,2)
      powER = (two*dnkk+one)/two
      do id = 1,nbas(isym)
      jd=id+nskip(isym)
      cl=two**(two*dnkk+onehalf)*(alpha(jd,isym)**power)/(ifact*spi)
      cs =dsqrt((two*dnkk-one)/(alpha(jd,isym)*
     :           (four*(dkap*dkap+dkap+dnkk)-one)))
      norm_large=dsqrt(cl/two)
      norm_small=norm_large*cs
      gl(1,jd)=0
      gs(1,jd)=0
      do ie = 2,n
      expon= dexp(-alpha(jd,isym)*r(ie)*r(ie))*rnk(ie)
      gl(ie,jd) = norm_large*expon
      gs(ie,jd) = norm_small*((dnkk+dkap)/r(ie)-
     :                   two*alpha(jd,isym)*r(ie))*expon
      qqf(ie,jd) = norm_large*((dnkk+dkap)/r(ie)-
     :                   two*alpha(jd,isym)*r(ie))*expon
c       print*,'gl,gs',gl(100,3),gs(100,3)
      enddo
      enddo
c       print*,'gl,gs',gl(100,3),gs(100,3)
      enddo

 200  format(/,2x,31hTOTAL NUMBER OF BASIS FUNCTION=i4,//)
 201  format(2x,'SYMMETRY   NO. OF BASIS FUNCTION   NSKIP',/)
 202  format(2x,i4,12x,i4,12x,i4)
      return
      end
c------------------------------------------------------------------
c     read input data
      subroutine readinp(mocc,method)
C------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c     include 'mpif.h'
      character*80 title
      character*14 SCFIN,OUTPUT
      character*2 atom,method
      common /datt/nsym,nocc,nbas(MNS),ind
     :  /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
      common /def1/amass,z
     :  /jvalue/orbj(MNOCC,MNS),iocc(MNBAS)
     :  /conv/crit,maxit,npower
     :  /jv/occ(MNOCC),orbe(MNB),iq(MNB)
     :  /print/jprint,iout
     :  /info/nocorb(MNOCC),nucorb(MNB)
     :  /nucop/nucopt
c    :  /UB/alpha1,beta1
c    :  /WT/alpha2,beta2,gamma,delta
     :  /ET1/alfa(MNS),beta(MNS),it_max
     :  /ET/alpha3(MNS),beta3(MNS)
     :  /freez/enmin(MNS),enmax(MNS)

     

      open(unit=IPUNCH,status='old',file='inputfile')
      read(IPUNCH,'(a)')SCFIN
      read(IPUNCH,'(a)')OUTPUT
      open (unit=STDIN,file=SCFIN,status='old')
      open (unit=STDOUT,file=OUTPUT,status='unknown')
      read(STDIN,'(a)')title
      read(STDIN,'(a)')atom
      read(STDIN,*)nsym,(nbas(i),i=1,nsym)
      read(STDIN,*)(nocorb(i),i=1,nsym)
      read(STDIN,*)(nak(i),i=1,nsym)
      read(STDIN,*)(occ(i),i=1,nsym)
      read(STDIN,'(a)')method
   
      alpha1=zero
      beta1=zero

      if(method.eq.'UB')then
      iopt=1
      endif
      if(method.eq.'WT')then
      iopt=2
      endif
      if(method.eq.'ET')then
      iopt=3
      endif
      print*,method
      goto  (1, 2, 3),iopt
   1   read(STDIN,*)alpha1,beta1
c  1   alpha1=alfa
c      beta1=beta
c      print*,'alpha1,beta1',alpha1,beta1

      goto 4
  2   read(STDIN,*)alpha2,beta2,gamma,delta 
      goto 4
c  3   read(STDIN,*)(alpha3(i),i=1,nsym)
c      read(STDIN,*)(beta3(i),i=1,nsym)
  3    do i=1,it_max
       alpha3(i)=alfa(i)
       beta3(i)=beta(i)
       enddo
c      print*,'alpha',(alpha3(i),i=1,it_max)
c      print*,'beta',(beta3(i),i=1,it_max)
  4   continue
      read(STDIN,*)icharg
      read(STDIN,*)maxit,npower
*  iout=0 is output without energy limit, otherwise with limit
      read(STDIN,*)jprint,iout
      read(STDIN,*)nucopt
      read(STDIN,*)amass,z

      if (iout.ne.0) then
      read(STDIN,*)(enmin(i),i=1,nsym)
      read(STDIN,*)(enmax(i),i=1,nsym)
      endif

      crit=(0.1d0)**npower
      iele=0
      mocc=0
      do isym=1,nsym
      nn=nocorb(isym)
      if(nn.ne.0)then
      do ibas=1,nocorb(isym)
      iele=iele+2*occ(isym)+1
      mocc=mocc+1
      enddo
      endif
      enddo
      call scfout
      write(STDOUT,9994)title
      write(STDOUT,301)
      call nucleus(atom,nnuc)
      write(STDOUT,302)nsym
      write(STDOUT,*)'BASIS   --> ',(nbas(i),i=1,nsym)
      write(STDOUT,*)'OCC ORB --> ',(nocorb(i),i=1,nsym)
      write(STDOUT,*)'KAPPA   --> ',(nak(i),i=1,nsym)
      write(STDOUT,*)'J VALUE --> ',(occ(i),i=1,nsym)
      write(STDOUT,9995)method
      goto (5,6,7),iopt
  5   write(STDOUT,310)alpha1,beta1
      go to 8
  6   write(STDOUT,311)alpha2,beta2,gamma,delta
      go to 8
  7   write(STDOUT,*)'ALPHA3  --> ',(alpha3(i),i=1,nsym)
      write(STDOUT,*)'BETA3   --> ',(beta3(i),i=1,nsym)
  8   continue
      write(STDOUT,303)icharg
      write(STDOUT,304)maxit, crit
      write(STDOUT,305)jprint
      write(STDOUT,306)nucopt
      write(STDOUT,309)
  301 format('>>>>>>>>>>>>>>> INPUT DATA >>>>>>>>>>>>>>>>>>>>>>>>')
  302 format('NUMBER OF SYMMETRY = ',i2)
  304 format('MAX ITER = ',i2,8x,'CONV. =',d10.4)
  303 format('CHARGE =',i2)
  305 format('JPRINT =',i2)
  306 format('NUCLEAR OPTIONAL = ',i2)
  309 format('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<',/)
  310 format('ALPHA1 = ',d10.4,8x,'BETA1 =',d10.4)
  311 format('ALPHA2 = ',d10.4,3x,'BETA2 =',d10.4,3x,'GAMMA = ',
     :   d10.4,3x,'DELTA =',d10.4)
9994  format(//,27x,a,//)
9995  format(//,25x,a,//)
      if(iele+icharg.ne.nnuc)then
      write(STDOUT,*)'>>>>>>>>> PROGRAMME TERMINATES DUE TO IO ERROR'
      write(STDOUT,*)' CHECK ATOMIC NUMBER, CHARGE AND OCCUPANCY'
      stop
      else
      z=dble(nnuc)
      endif
      open (unit=IWFN,file='wfn.dat',form='unformatted',
     :          status='unknown')
      open (unit=IWFN,file='wfn.dat',status='old')

      return
      end
c---------------------------------------------------------------
      subroutine nucparm
c---------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c     include 'mpif.h'
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
      common/cons/zero,half,tenth,one,two,three,ten
     :  /npar/parm(2),nparm
     :  /def1/amass,z
      common/wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
      common/quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)


      nparm=2
      pi = 3.141592653589790d 000
      tfm=2.3d0
      ainfcm = 0.529177249d-08
      fmtoau = 1.0d-13/ainfcm
      t=tfm*fmtoau
      const= 4.0d0* log (3.0d0)
      parm(2)=t/const
      rrmsfm=0.836d0*amass**(1.d0/3.d0)+0.57d0
      rrms = rrmsfm*fmtoau
      factor = rrms**2-(7.0d 00/5.0d 00)
     :                    *(pi**2)*(parm(2)**2)
      parm(1) =   sqrt (5.0d 00/3.0d 00)* dsqrt (factor)
      return
      end
c-----------------------------------------------------------
      subroutine findk(orba,orbb,orbp,orbq,ia,ib,ip,iq,
     $jeven,kkk)
c     this subroutine finds the common k value
c-----------------------------------------------------------
      implicit real*8 (a-h,o-z)
c     include 'mpif.h'
      real*8 jj
      dimension ja(MAXVAL),jb(MAXVAL),jp(MAXVAL),jq(MAXVAL),
     :          jeven(MAXVAL)

      do i=1,MAXVAL
          jeven(i)=0
      enddo

      iab=ia*ib
      ipq=ip*iq

      jminab=idint(dmax1(orba,orbb)-dmin1(orba,orbb))
      jminpq=idint(dmax1(orbp,orbq)-dmin1(orbp,orbq))

      jmaxab=idint(orba+orbb)
      jmaxpq=idint(orbp+orbq)

      ll=0
      kk=0
      do i=jminab,jmaxab
          jtot=jmaxab+i
          jj=float(jtot)-2.0d0*(float(jtot/2))
          if(jj.ne.0.and.iab.gt.0)then
          ll=ll+1
          ja(ll)=i
      endif
      if(jj.eq.0.and.iab.lt.0)then
      kk=kk+1
      jb(kk)=i
      endif
      enddo
      mm=0
      nn=0
      do i=jminpq,jmaxpq
      jtot=jmaxpq+i
      jj=float(jtot)-2.0*(float(jtot/2))
      if(jj.ne.0.and.ipq.gt.0)then
      mm=mm+1
      jp(mm)=i
      endif
      if(jj.eq.0.and.ipq.lt.0)then
      nn=nn+1
      jq(nn)=i
      endif
      enddo
      kkk=0
      do i=1,ll
      do j=1,mm
      if(ja(i).eq.jp(j))then
      kkk=kkk+1
      jeven(kkk)=ja(i)
      endif
      enddo
      enddo

      do i=1,ll
      do j=1,nn
      if(ja(i).eq.jq(j))then
      kkk=kkk+1
      jeven(kkk)=ja(i)
      endif
      enddo
      enddo
      do i=1,kk
      do j=1,mm
      if(jb(i).eq.jp(j))then
      kkk=kkk+1
      jeven(kkk)=jb(i)
      endif
      enddo
      enddo
      do i=1,kk
      do j=1,nn
      if(jb(i).eq.jq(j))then
      kkk=kkk+1
      jeven(kkk)=jb(i)
      endif
      enddo
      enddo
      return
      end

      subroutine integral(nbasis)
c-------------------------------------------------------------------
c This evaluates the direct and exchange matrix elements numericaly.
c It computes submatrices skll,skss and c skls. the submatrix skls 
c is not computed as the fock matrix is symmetric. P(i,j) is the
c density matrix.
c-------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c     include 'mpif.h'
      dimension momj(MAXVAL),jq(MNS)
      common/cons/zero,half,tenth,one,two,three,ten
     : /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     : /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     : /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     : /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     :  /diag/eigv(MNB,MNB),eng(MNB)
     : /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
     : /double/df_two(MNB,MNB)
     : /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     : /datt/nsym,nocc,nbas(MNS),ind
     :  /jvalue/orbj(MNOCC,MNS),iocc(MNBAS)
     :  /info/nocorb(MNOCC),nucorb(MNB)
     :  /skip/nskip(MNS)
     :  /dens/p(MNB,MNB)
     :  /kapp/kappa(MNB),kappc(MNB),iqc(MNOCC)
     :  /jv/occ(MNOCC),orbe(MNB),iq(MNB)
     :  /info1/occk(MNOCC,MNS)
       data jq/1,-1,1,-1,1,-1,1,-1,1,-1,1/
      do ia = 1,MNB
      do ib = 1,MNB
      df_two(ib,ia) = zero
      enddo
      enddo

      mtp=n
      do isym=1,nsym
      do ia=1,nbas(isym)
      ka=ia+nskip(isym)
      ja=ia+2*nskip(isym)
      index1=ia+nbas(isym)+2*nskip(isym)
      do ib=1,ia
      index2=ib+nbas(isym)+2*nskip(isym)
      jb=ib+2*nskip(isym)
      kb=ib+nskip(isym)

       df_two(ja,jb)=zero
       df_two(index1,index2)=zero
       df_two(ja,index2)=zero
       df_two(index1,jb)=zero
      
      do jsym=1,nsym
      nn=nocorb(jsym)
      if(nn.ne.0)then
      do ic=1,nocorb(jsym)
      orbc=occ(jsym)
      orba=orbe(ka)
      orbb=orbe(kb)
      iaa=iq(ka)
      ibb=iq(kb)
      icc=jq(jsym)
      kc=ic+nskip(jsym)
c
c     direct part
c
c     find k
      call findk(orba,orbb,orbc,orbc,iaa,ibb,icc,icc,momj,lmax)
      if(lmax.ne.0) then
      do  lk=1,lmax
      ll=momj(lk)
      if(ll.eq.0)then
      call slater1(ka,kc,kb,kc,ll,rkll,rkss)


      df_two(ja,jb)=df_two(ja,jb)+rkll*(2.d0*orbc+1.d0)

      df_two(jb,ja)=df_two(ja,jb)

      df_two(index1,index2)=df_two(index1,index2)+rkss*(2.d0*orbc+1.d0)
      
      df_two(index2,index1)=df_two(index1,index2)
      endif
      enddo
      endif
c
c     exchange part
c
c     find k
      call findk(orba,orbc,orbc,orbb,iaa,icc,icc,ibb,momj,kmax)
      if(kmax.ne.0)then
      do ik=1,kmax
      kk=momj(ik)
      fact=(2.d0*orbc+1.d0)*clrx (kappa(ka),kk,kappa(kc))**2

      call slater2(ka,kc,kc,kb,kk,rkll,rkss,rkls,rksl)
      df_two(ja,jb)=df_two(ja,jb)-rkll*fact
      df_two(jb,ja)=df_two(ja,jb)
      df_two(index1,index2)=df_two(index1,index2)-rkss*fact
      df_two(index2,index1)=df_two(index1,index2)
      df_two(ja,index2)=df_two(ja,index2) -rkls*fact
      df_two(index2,ja)=df_two(ja,index2)
      df_two(index1,jb)=df_two(index1,jb) -rksl*fact
      df_two(jb,index1)=df_two(index1,jb)
      enddo
      endif
      enddo
      endif
      enddo

      enddo
      enddo
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine nucleus(tex,nmas)
      parameter (kerne = 105)
      character*2 text,tex
      dimension text(kerne)

      data text/'  ',' H','HE','LI','BE',' B',' C',' N',' O',' F','NE',
     &'NA','MG','AL','SI',' P',' S','CL','AR',' K','CA','SC','TI',' V',
     &'CR','MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR',
     &'RB','SR',' Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN',
     &'SN','SB','TE',' I','XE','CS','BA','LA','CE','PR','ND','PM','SM',
     &'EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA',' W','RE',
     &'OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN','FR','RA',
     &'AC','TH','PA',' U','NP','PU','AM','CM','BK','CF','ES','FM','MD',
     &'NO','LR','04'/

      do 1 nl=1,kerne
      if  (tex.eq.text(nl)) goto 2
    1 continue
      write(STDOUT,201) tex,nmas
  201 format(//,' ATOM ',a2,i4,'  UNKNOWN')
      stop
    2 nmas=nl-1
      write(STDOUT,202) tex,nmas
  202 format(//,'    ATOM =   ',a2,8x,' ATOMIC NUMBER = ',i3)
      return
      end
c------------------------------------------------------------------
      subroutine twoe(isym,nbasisi,iter)
c-------------------------------------------------------------------
c This evaluates the direct and exchange matrix elements numericaly.
c It computes submatrices skll,skss and c skls. the submatrix skls 
c is not computed as the fock matrix is symmetric. P(i,j) is the
c density matrix.
c------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c     include 'mpif.h'
c     real*8 df_two_local
      dimension momj(MAXVAL),jq(MNS)
c    :,df_two_local(MNB,MNB)
      common/cons/zero,half,tenth,one,two,three,ten
     : /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     : /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     : /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     : /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
     :  /diag/eigv(MNB,MNB),eng(MNB)
     : /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
     : /double/df_two(MNB,MNB)
     : /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     : /datt/nsym,nocc,nbas(MNS),ind
     :  /jvalue/orbj(MNOCC,MNS),iocc(MNBAS)
     :  /info/nocorb(MNOCC),nucorb(MNB)
     :  /skip/nskip(MNS)
     :  /dens/p(MNB,MNB)
     :  /kapp/kappa(MNB),kappc(MNB),iqc(MNOCC)
     :  /jv/occ(MNOCC),orbe(MNB),iq(MNB)
     :  /info1/occk(MNOCC,MNS)
       data jq/1,-1,1,-1,1,-1,1,-1,1,-1,1/
c     call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
c     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      do ia = 1,MNB
      do ib = 1,MNB
      df_two(ib,ia) = zero
c     df_two_local(ib,ia) = zero
      enddo
      enddo
      fact=0.0
c      print*,'df_two',df_two(10,10)
      mtp=n
      do ia=1,nbas(isym)
      ka=ia+nskip(isym)
      ja=ia+2*nskip(isym)
      index1=ia+nbas(isym)+2*nskip(isym)
      do ib=1,ia
      index2=ib+nbas(isym)+2*nskip(isym)
      jb=ib+2*nskip(isym)
      kb=ib+nskip(isym)

      df_two(ja,jb)=zero
      df_two(index1,index2)=zero
      df_two(ja,index2)=zero
      df_two(index1,jb)=zero

      do jsym=1,nsym
      nn=nocorb(jsym)
      if(nn.ne.0)then
      do ic=1,nocorb(jsym)
      orbc=occ(jsym)
      orba=orbe(ka)
      orbb=orbe(kb)
      iaa=iq(ka)
      ibb=iq(kb)
      icc=jq(jsym)
      kc=ic+nskip(jsym)
c
c     direct part
c
c     find k
      call findk(orba,orbb,orbc,orbc,iaa,ibb,icc,icc,momj,lmax)
      if(lmax.ne.0) then
      do  lk=1,lmax
      ll=momj(lk)
      if(ll.eq.0)then
      call slater1(ka,kc,kb,kc,ll,rkll,rkss)
c     print*,'rkll,rkss',rkll,rkss 
      df_two(ja,jb)=df_two(ja,jb)+rkll*(2.d0*orbc+1.d0)

      df_two(jb,ja)=df_two(ja,jb)

      df_two(index1,index2)=df_two(index1,index2)
     :  +rkss*(2.d0*orbc+1.d0)
      
      df_two(index2,index1)=df_two(index1,index2)
 
c     if (iter.eq.2.and.ja.eq.3.and.jb.eq.3)
c    : print*,'dfff',df_two(ja,jb)

      endif
      enddo
      endif
c
c     exchange part
c
c     find k
      call findk(orba,orbc,orbc,orbb,iaa,icc,icc,ibb,momj,kmax)
      if(kmax.ne.0)then
      do ik=1,kmax
      kk=momj(ik)
      fact=(2.d0*orbc+1.d0)*clrx (kappa(ka),kk,kappa(kc))**2

      call slater2(ka,kc,kc,kb,kk,rkll,rkss,rkls,rksl)
      df_two(ja,jb)=df_two(ja,jb)-rkll*fact
      df_two(jb,ja)=df_two(ja,jb)
      df_two(index1,index2)=df_two(index1,index2)
     :  -rkss*fact
      df_two(index2,index1)=df_two(index1,index2)
      df_two(ja,index2)=df_two(ja,index2) -rkls*fact
      df_two(index2,ja)=df_two(ja,index2)
      df_two(index1,jb)=df_two(index1,jb) -rksl*fact
      df_two(jb,index1)=df_two(index1,jb)
      enddo
      endif
      enddo
      endif
      enddo

      enddo
      enddo
c      print*,'df_two',df_two(10,10)
c            call MPI_ALLREDUCE(df_two_local,df_two,MNB*MNB,MPI_REAL8
c    : ,MPI_SUM,MPI_COMM_WORLD,ierr)

c     do i= 1,MNB
c      do j=1,MNB
c      df_two_local(i,j)=0.0
c      enddo
c     enddo

      return
      end
c***************************************************************
      subroutine eigen (a,b,nsub,valu,msub,anorm,nnmax)
      parameter(nmax=1600)
c     eigenvalues and eigenvectors of a real symmetric matrix
c     a ----- secular determinant (input)
c             eigenvetors  column-wise (output)
c             a is destroyed
c     b ----- nsub x nsub array for intermidiate storage
c     nsub ----- order of a (input)
c                must be equal to or greater than 1
c     valu ----- eigenvalues in decreasing order (output)
c     msub -----  number of eigenvalues (input)
c                 if equal to -1, eigenvectors are not computed
c     anorm ----- matrix norm which is a sqrt of sum of a(i,j)**2
c                                        over all i & j (output)
c     nmax ----- row dimension of a & b
c
      implicit real*8 (a-h,o-z)
c     include 'mpif.h'
      integer index
      dimension a(nnmax,nnmax), b(nnmax,nnmax), valu(nnmax)
      dimension diag(nmax), superd(nmax), wvec(nmax),pvec(nmax),
     :          qvec(nmax), vall(nmax), q(nmax), u(nmax),
     :          index(nmax), factor(nmax), v(nmax), t(nmax,3)
      common /f202/ diag, superd, wvec, pvec, t
      equivalence (wvec(1),vall(1),factor(1),u(1)),
     1            (pvec(1),qvec(1),q(1),v(1))
      data e1/2.3d-15/
c     householder similarity transformation to co-diagonal form
c     ----------- ---------- -------------- -- ----------- ----
      lpout = 0
      anorm=0
      n=nsub
      m=msub
      do 10 i = 1,n
      valu(i) =0.0d+00
   10 continue
      if (n.lt.0) go to 220
      if (m.le.0.or.n.le.2) go to 50
      do 40 i=2,n
      do 40 j=2,n
      b(j,i)=0.0d0
      if (i.eq.j) b(j,i)=1.0d0
   40 continue
   50 do 200 i=1,n
      i1=i+1
      i2=i1+1
      if (i2.gt.n) go to 160
      sum=0.0d0
      do 70 j=i2,n
      sum=sum+a(j,i)**2
      wvec(j)=0.0
   70 continue
      if (sum.eq.0.0d0) go to 160
      j=i1
      temp=a(j,i)
      sum=dsqrt(sum+ temp  **2)
      a(j,i)=-dsign(sum, temp  )
      wvec(j )=dsqrt(     1.0d0+dabs( temp  )/sum)
      div=dsign(    wvec(j )*sum, temp  )
      do 85 j=i2,n
   85 wvec(j)=a(j,i)/div
      scalar=0.0d0
      do 95 j=i1,n
      pvec(j)=0.0d0
      qvec(j)=0.0d0
      do 90 k=i1,n
   90 pvec(j)=pvec(j)+a(k,j)*wvec(k)
      scalar=scalar+pvec(j)*wvec(j)
   95 continue
      scalar=scalar/2.0d0
      do 120 j=i1,n
      qvec(j)=pvec(j)-scalar*wvec(j)
      do 120 k=i1,j
      a(k,j)=a(k,j)-(wvec(k)*qvec(j)+wvec(j)*qvec(k))
      a(j,k)=a(k,j)
  120 continue
      if (m.le.0) go to 160
      do 150 k=2,n
      temp=0.0d0
      do 140 j=i1,n
  140 temp=temp+wvec(j)*b(j,k)
      do 150 j=i1,n
      b(j,k)=b(j,k)-wvec(j)*temp
  150 continue
  160 j=i
      diag(i)=a(j,i)
      if (i.ne.n) superd(i)=a(j+1,i)
  200 continue
c     givens eigenvalue iteration from sturm chain of co-diagonal minors
  220 n=iabs (n)
      m=iabs (m)
c     calculate norm of matrix and initialize eigenvalue bounds
      anorm2=0.0
      anorm2=diag(1)**2
      if (n.eq.1) go to 235
      do 230 l=2,n
      q(l-1)=superd(l-1)**2
      anorm2=diag(l)**2+q(l-1)+q(l-1)+anorm2
  230 continue
  235 anorm=0.0
      anorm=dsqrt(anorm2)
      if (m.eq.0) return
      do 240 l=1,m
      valu(l)=anorm
      vall(l)=-anorm
  240 continue
      eps1=anorm*e1
      if (eps1.eq.0.0d0) return
  250 do 570 l=1,m
  260 tau=(valu(l)+vall(l))/2.0d0
      if (    (tau-vall(l)).le.e1*(dabs(tau)+e1)) go to 570
      match=0
      t2=0.0d0
      t1=1.0d0
      do 450 l1=1,n
      p=diag(l1)-tau
      if (t2.ne.0.0d0) go to 330
  300 t1=dsign(1.0d0,t1)
  330 if (t1.ne.0.0d0) go to 400
      t0=-dsign(1.0d0,t2)
      t2=0.0d0
      if (q(l1-1).ne.0.0d0) go to 410
      go to 300
  400 t0=p-q(l1-1)*t2/t1
      t2=1.0d0
  410 if (t0) 440, 420, 430
  420 t2=t1
      if (t2) 440, 430, 430
  430 match=match+1
  440 t1=t0
  450 continue
      do 530 l1=l,m
      if (l1.le.match) go to 500
      if (valu(l1).le.tau) go to 260
      valu(l1)=tau
      go to 530
  500 vall(l1)=tau
  530 continue
      go to 260
  570 continue
c     eigenvectors of co-diagonal symmetric matrix--inverse iteration
      m=msub
      if (m.lt.0) return
      do 970 i=1,m
      if (i.eq.1) go to 725
      if (valu(i-1)-valu(i).lt.1.0e+4*eps1) go to 730
  725 i1=-1
  730 i1=i1+1
      do 760 l=1,n
      v(l)=eps1*(i1*l+1)
      t(l,2)=diag(l)-valu(i)
      if (l.lt.n) go to 740
      t(l,3)=0.0d0
      go to 760
  740 t(l,3)=superd(l)
      if (t(l,3).eq.0.0d0) t(l,3)=eps1
      t(l+1,1)=t(l,3)
  760 continue
      do 820 j=1,n
      t(j,1)=t(j,2)
      t(j,2)=t(j,3)
      t(j,3)=0.0d0
      go to 780
  770 t(j,1)=eps1
c**** replace function name  in next      statement  in single precision
  780 vtemp=dabs(t(j,1))
      if (vtemp.lt.eps1) go to 770
      if (j.eq.n) go to 820
      index(j)=0
      if (dabs(t(j+1,1)).le.vtemp) go to 810
      index(j)=1
      do 800 k=1,3
      vtemp=t(j,k)
      t(j,k)=t(j+1,k)
      t(j+1,k)=vtemp
  800 continue
  810 vtemp   =t(j+1,1)/t(j,1)
      factor(j)=vtemp
      t(j+1,2)=t(j+1,2)-  vtemp  *t(j,2)
      t(j+1,3)=t(j+1,3)-  vtemp  *t(j,3)
  820 continue
      iter=1
      if (i1.gt.0) go to 920
  860 do 870 l1=1,n
      l=n+1-l1
      v( l )=(v( l )-t(l,2)*v( l+1 )-t(l,3)*v( l+2 ))/t(l,1)
  870 continue
      go to (875,920), iter
  875 iter=2
      if (n.eq.1) go to 860
  880 do 910 l=2,n
      if (index(l-1).eq.0) go to 900
      vtemp=v(l-1)
      v(l-1)=v(l)
      v(l)=vtemp
  900 v(l)=v(l)-factor(l-1)*v(l-1)
  910 continue
      go to 860
  920 if (i1.eq.0) go to 945
      do 940 l1=1,i1
      k=i-l1
      vtemp=0.0d0
      do 930 j=1,n
  930 vtemp=vtemp+a(j,k)*v(j)
      do 940 j=1,n
  940 v(j)=v(j)-a(j,k)*vtemp
      go to (880,945), iter
  945 vnorm2=0.0d0
      do 950 l=1,n
  950 vnorm2=vnorm2+v( l )**2
      vnorm=dsqrt(vnorm2)
      do 960 j=1,n
  960 a(j,i)=v( j )/vnorm
  970 continue
      n=nsub
      if (n.le.2) return
      do 990 i=1,m
      do 980 k=2,n
      u(k)=0.0d0
      do 980 j=2,n
  980 u(k)=u(k)+b(j,k)*a(j,i)
      do 990 j=2,n
  990 a(j,i)=u(j)
      return
      end
c***********************************************************************
      subroutine slater1(ia,ib,ic,id,k,rkll,rkss)

*   the value of this  function is the  slater integral  as normally   *
*   defined in terms of the four sets of quantum numbers  a,b,c,d.     *
*                                                                      *
*   subroutines called: quad, yzk.                                     *
*                                                                      *
c***********************************************************************
*
      implicit real*8(a-h, o-z)
c     include 'mpif.h'
c     dimension rhop(MN),rttk(MN),wk(MN),temp(MN),zk(MN)
      common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)

      mtp=n
      do i=1,mtp
         tmp1(i)=0
         tmp3(i)=0
c       rkll=00.0
c       rkss=0.0
      enddo

      call yzk1(k,ib,id)
      do i = 2,mtp
         tmp3(i)=gl(i,ia)*gl(i,ic)*rpor(i)*tmp1(i)
      enddo
      call quad (tmp3,mtp,rkll)

      do i = 2,mtp
         tmp3(i)=gs(i,ia)*gs(i,ic)*rpor(i)*tmp1(i)
      enddo

      call quad (tmp3,mtp,rkss)

      return
      end
c***********************************************************************
      subroutine slater2(ia,ib,ic,id,k,rkll,rkss,rkls,rksl)

*   the value of this  function is the  slater integral  as normally   *
*   defined in terms of the four sets of quantum numbers  a,b,c,d.     *
*                                                                      *
*   subroutines called: quad, yzk.                                     *
*                                                                      *
c***********************************************************************
*
      implicit real*8(a-h, o-z)
c     include 'mpif.h'
c     dimension rhop(MN),rttk(MN),wk(MN),temp(MN),zk(MN)
      common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)

      mtp=n
      do i=1,mtp
         tmp1(i)=0
         tmp2(i)=0
         tmp3(i)=0
         tmp4(i)=0
      enddo

      call yzk2(k,ib,id)
      do i = 2,mtp
         tmp3(i)=gl(i,ia)*pf(i,ic)*rpor(i)*tmp1(i)
c      if (i.eq.100.and.ia.eq.3) print*,'pf',pf(i,ia),gl(i,ia)
      enddo
      call quad (tmp3,mtp,rkll)

      do i = 2,mtp
         tmp3(i)=gs(i,ia)*qf(i,ic)*rpor(i)*tmp2(i)
      enddo

      call quad (tmp3,mtp,rkss)
 
      do i = 2,mtp
         tmp3(i)=gs(i,ia)*qf(i,ic)*rpor(i)*tmp1(i)
      enddo

      call quad (tmp3,mtp,rksl)

      do i = 2,mtp
         tmp4(i)=gl(i,ia)*pf(i,ic)*rpor(i)*tmp2(i)
      enddo

      call quad (tmp4,mtp,rkls)

      return
      end
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine yzk1(k,i,j)
c-----------------------------------------------------------------------
c   this subroutine evaluates hartree y- and z-functions:              *
c                                                                      *
c               (k)            (k)           (k)                       *
c              y   (i,j;r) =  z   (i,j;r) + w   (i,j;r)                *
c                                                                      *
c   where                                                              *
c                                                                      *
c    (k)                                                               *
c   z   (i,j;r) =  i ( (s/r)   (p (s)*p (s) + q (s)*q (s)) ; 0 - r )   *
c                                i     j       i     j                 *
c                                                                      *
c   where                                                              *
c                                                                      *
c    (k)                    k+1                                        *
c   w   (i,j;r) =  i ( (r/s)   (p (s)*p (s) + q (s)*q (s)) ; r -       *
c                                i     j       i     j    infinity )   *
c                                                                      *
c   where  i ( g(r,s) ; range )  denotes the integral of g(r,s) over   *
c   range  in  s .  the y-function is tabulated in  common/tatb/  in   *
c   array  tb , the z-function in array ta .                           *
c-----------------------------------------------------------------------
      implicit real*8(a-h, o-z)
c     include 'mpif.h'
      dimension rhop(MN),rttk(MN),wk(MN),temp(MN),zk(MN)
      common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
      accy=1.0d-12
c----------------------------------------------------------------------
c   for  k > 0  compute  r^k   and store in  rttk
c----------------------------------------------------------------------
      if (k .gt. 0) then
         do 1 ii = 2,n
            rttk(ii) = r(ii)**k
    1    continue
      endif
c---------------------------------------------------------------------
c     determine maximum tabulation point as location beyond which rhop (see 
c comment statements below) would be zero; determine other important locations
c--------------------------------------------------------------------
      mtp = min(mf(i),mf(j))
      mtpp1 = mtp+1
      mtpp3 = mtp+3
      mtpp4 = mtp+4
c--------------------------------------------------------------------
c   compute rp(s)*(p_i (s)*p_j (s)+q_i (s)*q_j (s)) and store in rhop
c--------------------------------------------------------------------
      do 2 ii = 2,mtp
         rhop(ii) = rp(ii)*(pf(ii,i)*pf(ii,j)+qf(ii,i)*qf(ii,j))
    2 continue
c--------------------------------------------------------------------
c                   fill array temp with r**k * rhop
c--------------------------------------------------------------------
      temp(1) = zero
      if (k .eq. 0) then
         do 3 ii = 2,mtp
            temp(ii) = rhop(ii)
    3    continue
      else
         do 4 ii = 2,mtp
            temp(ii) = rttk(ii)*rhop(ii)
    4    continue
      endif
c---------------------------------------------------------------------
c                 set an additional four points to zero
c---------------------------------------------------------------------
      do 5 ii = mtpp1,mtpp4
         temp(ii) = zero
    5 continue
c-----------------------------------------------------------------------
c           compute the first few values of  r^k  * zk  using semi-open
c                            newton-cotes formulae
c-----------------------------------------------------------------------
      zk(1) = zero
      do 7 ii = 2,4
         sum = zero
         do 6 kk = 2,5
            sum = sum+cnc5c(kk,ii)*temp(kk)
    6    continue
         zk(ii) = sum
    7 continue
c-------------------------------------------------------------------
c   compute remainder of r^k  * zk: march out to mtp+3; use closed
c                       newton-cotes formula
c-------------------------------------------------------------------
      do 8 ii = 5,mtpp3
         zk(ii) = zk(ii-4)+c1*(temp(ii-4)+temp(ii  ))
     :                    +c2*(temp(ii-3)+temp(ii-1))
     :                    +c3* temp(ii-2)
    8 continue
c------------------------------------------------------------------------------
c determine the asymptotic value of  r^k * z^(k) correction to z^(0): in the 
c      manner of  c froese fischer,the hartree-fock method for atoms, 
c                  john wiley & sons,new york, 1977, p 235.
c------------------------------------------------------------------------------
      if (k .eq. 0) then
         if (i .eq. j) then
            zklim = one
         else
            zklim = zero
         endif
         do 10 kk = mtpp3,mtp,-1
            dif = zk(kk)-zklim
            if (abs (dif) .gt. accy) then
               do 9 ii = kk,2,-4
                  zk(ii) = zk(ii)-dif
    9          continue
            endif
   10    continue
      else
         zklim = zk(mtpp3)
      endif
c------------------------------------------------------------------------------
c                tabulate  zk  for entire internal lgrid
c------------------------------------------------------------------------------
      if (k .eq. 0) then
         do 11 ii = mtpp4,n
            zk(ii) = zklim
   11    continue
      else
         do 12 ii = 2,mtpp3
            zk(ii) = zk(ii)/rttk(ii)
   12    continue
         do 13 ii = mtpp4,n
            zk(ii) = zklim/rttk(ii)
   13    continue
      endif
c------------------------------------------------------------------------------
c                       start array wk / r**(k+1)
c------------------------------------------------------------------------------
      np4 = n+4
      do 14 ii = np4,mtpp1,-1
         wk(ii) = zero
   14 continue
c------------------------------------------------------------------------------
c                   compute  r^(k+1) and store in rttk
c------------------------------------------------------------------------------
      if (k .gt. 0) then
         do 15 ii = 2,n
            rttk(ii) = rttk(ii)*r(ii)
   15    continue
      endif
c------------------------------------------------------------------------------
c        fill array temp with rhop / r**(k+1) ; set temp(1) = zero
c                         to avoid 0/0 case
c------------------------------------------------------------------------------
      temp(1) = zero
      if (k .eq. 0) then
         do 16 ii = 2,mtp
            temp(ii) = rhop(ii)/r(ii)
   16    continue
      else
         do 17 ii = 2,mtp
            temp(ii) = rhop(ii)/rttk(ii)
   17    continue
      endif
c------------------------------------------------------------------------------
c         compute remainder of wk / r**(k+1): march in to the origin
c------------------------------------------------------------------------------
      do 18 ii = mtp,2,-1
         wk(ii) = wk(ii+4)+c1*(temp(ii  )+temp(ii+4))
     :                    +c2*(temp(ii+1)+temp(ii+3))
     :                    +c3*(temp(ii+2))
   18 continue
      wk(1) = zero
c------------------------------------------------------------------------------
c                               compute wk
c------------------------------------------------------------------------------
      if (k .eq. 0) then
         do 19 ii = 2,mtp
            wk(ii) = wk(ii)*r(ii)
   19    continue
      else
         do 20 ii = 2,mtp
            wk(ii) = wk(ii)*rttk(ii)
   20    continue
      endif
c------------------------------------------------------------------------------
c                           assemble solution
c------------------------------------------------------------------------------
      tmp1(1) = zero
      do 21 ii = 2,n
         tmp1(ii) = zk(ii)+wk(ii)
   21 continue
*
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine yzk2(k,i,j)
c------------------------------------------------------------------------------
c   this subroutine evaluates hartree y- and z-functions:              *
c               (k)            (k)           (k)                       *
c              y#   (i,j;r) =  z#   (i,j;r) + w#   (i,j;r)             *
c                                                                      * 
c   where # can be p or q and                                          *
c    (k)                                                               *
c   z#   (i,j;r) =  i ( (s/r)   # (s)*# (s) ; 0 - r )                  *
c                                i     j                               *
c    (k)                    k+1                                        *
c   w#   (i,j;r) =  i ( (r/s)    # (s)*# (s) ; r - infinity)           *
c                                 i     j                              *
c                                                                      *
c   where  i ( g(r,s) ; range )  denotes the integral of g(r,s) over   *
c   range  in  s .  the yp-function is tabulated in  common/tatb/  in  *
c   array  ta , the yq-function in array tb . for the direct part the  *
c   two araays yp and yq are summed and tabulated in array ta.         *
c------------------------------------------------------------------------------
      implicit real*8(a-h, o-z)
c     include 'mpif.h'
      dimension rhopp(MN),wkp(MN+5),tempp(MN),zkp(MN),
     :          rhopq(MN),wkq(MN+5),tempq(MN),zkq(MN),rttk(MN)
      common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
c------------------------------------------------------------------------------
c           for  k > 0  compute  r^k   and store in  rttk
c------------------------------------------------------------------------------
      if (k .gt. 0) then
        do 1 ii = 2,n
          rttk(ii) = r(ii)**k
    1   continue
      endif
c------------------------------------------------------------------------------
c   determine maximum tabulation point as location beyond which rhop  (see 
c     comment statements below) would be zero; determine other important 
c                               locations
c------------------------------------------------------------------------------
      mtp = mf(i)
      mtpp1 = mtp+1
      mtpp3 = mtp+3
      mtpp4 = mtp+4
c------------------------------------------------------------------------------
c   compute rp(s)* p_i (s)*gl_j (s) and store in rhopp
c                  q_i     gs_j
c------------------------------------------------------------------------------
      do 2 ii = 2,mtp
         rhopp(ii) = rp(ii)*pf(ii,i)*gl(ii,j)
         rhopq(ii) = rp(ii)*qf(ii,i)*gs(ii,j)
    2 continue
c------------------------------------------------------------------------------
c                   fill array tempp with r**k * rhopp
c                              tempq             rhopq
c------------------------------------------------------------------------------
      tempp(1) = zero
      tempq(1) = zero
      if (k .eq. 0) then
        do 3 ii = 2,mtp
          tempp(ii) = rhopp(ii)
          tempq(ii) = rhopq(ii)
    3   continue
      else
        do 4 ii = 2,mtp
          tempp(ii) = rttk(ii)*rhopp(ii)
          tempq(ii) = rttk(ii)*rhopq(ii)
    4   continue
      endif
c------------------------------------------------------------------------------
c                   set an additional four points to zero
c------------------------------------------------------------------------------
      do 5 ii = mtpp1,mtpp4
        tempp(ii) = zero
        tempq(ii) = zero
    5 continue
c------------------------------------------------------------------------------
c compute first few values of  r^k*zk using semi-open newton-cotes formulae
c------------------------------------------------------------------------------
      zkp(1) = zero
      zkq(1) = zero
      do 7 ii = 2,4
        sump = zero
        sumq = zero
        do 6 kk = 2,5
          sump = sump+cnc5c(kk,ii)*tempp(kk)
          sumq = sumq+cnc5c(kk,ii)*tempq(kk)
    6   continue
        zkp(ii) = sump
        zkq(ii) = sumq
    7 continue
c------------------------------------------------------------------------------
c        compute remainder of r^k*zk: march out to mtp+3; use closed
c                         newton-cotes formula
c------------------------------------------------------------------------------
      do 8 ii = 5,mtpp3
         zkp(ii) = zkp(ii-4)+c1*(tempp(ii-4)+tempp(ii  ))
     :                      +c2*(tempp(ii-3)+tempp(ii-1))
     :                      +c3* tempp(ii-2)
         zkq(ii) = zkq(ii-4)+c1*(tempq(ii-4)+tempq(ii  ))
     :                      +c2*(tempq(ii-3)+tempq(ii-1))
     :                      +c3* tempq(ii-2)
    8 continue
c------------------------------------------------------------------------------
c   determine the asymptotic value of  r^k*z^(k) correction to  z^(k)   
c         : in the manner of  c froese fischer,the hartree-fock method for 
c            atoms, john wiley & sons, new york, 1977, p 235.
c------------------------------------------------------------------------------
      if (k .eq. 0) then
        zklimp = zero
        zklimq = zero
      else
        zklimp = zkp(mtpp3)
        zklimq = zkq(mtpp3)
      endif
c------------------------------------------------------------------------------
c                     tabulate  zk  for entire internal lgrid
c------------------------------------------------------------------------------
      if (k .eq. 0) then
        do 12 ii = mtpp4,n
          zkp(ii) = zklimp
          zkq(ii) = zklimq
   12   continue
      else
        do 13 ii = 2,mtpp3
          zkp(ii) = zkp(ii)/rttk(ii)
          zkq(ii) = zkq(ii)/rttk(ii)
   13   continue
        do 14 ii = mtpp4,n
          zkp(ii) = zklimp/rttk(ii)
          zkq(ii) = zklimq/rttk(ii)
   14   continue
      endif
c------------------------------------------------------------------------------
c                          start array wk / r**(k+1)
c------------------------------------------------------------------------------
         np4 = n+4
         do 15 ii = np4,mtpp1,-1
            wkp(ii) = zero
            wkq(ii) = zero
   15    continue
c------------------------------------------------------------------------------
c               compute  r^(k+1)  and store in rttk
c------------------------------------------------------------------------------
         if (k .gt. 0) then
         do 16 ii = 2,n
            rttk(ii) = rttk(ii)*r(ii)
   16    continue
         endif
c------------------------------------------------------------------------------
c        fill array temp with rhop / r**(k+1) ; set temp(1) = zero
c                          to avoid 0/0 case
c------------------------------------------------------------------------------
         tempp(1) = zero
         tempq(1) = zero
         if (k .eq. 0) then
         do 17 ii = 2,mtp
            tempp(ii) = rhopp(ii)/r(ii)             
            tempq(ii) = rhopq(ii)/r(ii)
   17    continue
         else
         do 18 ii = 2,mtp
            tempp(ii) = rhopp(ii)/rttk(ii)
            tempq(ii) = rhopq(ii)/rttk(ii)
   18    continue
         endif
c------------------------------------------------------------------------------
c         compute remainder of wk / r**(k+1): march in to the origin
c------------------------------------------------------------------------------
      do 19 ii = mtp,2,-1
         wkp(ii) = wkp(ii+4)+c1*(tempp(ii  )+tempp(ii+4))
     :                    +c2*(tempp(ii+1)+tempp(ii+3))
     :                    +c3*(tempp(ii+2))
         wkq(ii) = wkq(ii+4)+c1*(tempq(ii  )+tempq(ii+4))
     :                    +c2*(tempq(ii+1)+tempq(ii+3))
     :                    +c3*(tempq(ii+2))
   19 continue
      wkp(1) = zero
      wkq(1) = zero
c------------------------------------------------------------------------------
c                                 compute wk
c------------------------------------------------------------------------------
      if (k .eq. 0) then
         do 20 ii = 2,mtp
            wkp(ii) = wkp(ii)*r(ii)
            wkq(ii) = wkq(ii)*r(ii)
   20    continue
      else
         do 21 ii = 2,mtp
            wkp(ii) = wkp(ii)*rttk(ii)
            wkq(ii) = wkq(ii)*rttk(ii)
   21    continue
      endif
c-----------------------------------------------------------
c                             assemble solution
c----------------------------------------------------------
      tmp1(1) = zero
      tmp2(1) = zero
      do 22 ii = 2,n
        tmp1(ii) = zkp(ii)+wkp(ii)
        tmp2(ii) = zkq(ii)+wkq(ii) 
   22 continue
      return
*
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine scfiter(nbasis,indt,isk)
c------------------------------------------------------------------
c     This performs the scf iteration.
c------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c     include 'mpif.h'
      logical set
      character*3 nh_i
      real*8 df_mat,swap
      dimension can(MNB,MNB)
      dimension  drs(MNB,MNB),sa(MNB),df_mat(MNB,MNB)
      common/basis/gl(MN,MNBAS),gs(MN,MNBAS),alpha(MNBAS,MNS),mg(MNBAS)
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /cons/zero,half,tenth,one,two,three,ten
     :  /overlap/sss(MNB,MNB),sll(MNB,MNB)
     :  /adat/occj(MNOCC),occl(MNOCC),unocl(MNS),unocj(MNS)
     :  /quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)
     :  /diag/eigv(MNB,MNB),eng(MNB)
     :  /datt/nsym,nocc,nbas(MNS),ind
     :  /single/df_single(MNB,MNB)
     :  /double/df_two(MNB,MNB)
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
     :  /conv/crit,maxit,npower
     :  /print/jprint,iout
     :  /jv/occ(MNOCC),orbe(MNB),iq(MNB)
     :  /skip/nskip(MNS)
     :  /dens/p(MNB,MNB)
     :  /ovlp/emt(MNB,MNB)
     :  /info1/occk(MNOCC,MNS)
     :  /info/nocorb(MNOCC),nucorb(MNB)
     :  /vector/vec(MNB,MNB),oe(MNB,MNS)
     :  /info3/df_t(MNB,MNB)
      common/optim/gauss1(200),gauss2(200),gauss3(200)
     :  /ET1/alfa(MNS),beta(MNS),it_max
     :  /totbasis/itot,norb(200),ida(13),icom(13),idat(13)

      dimension a(MNB,MNB),st(MNB,MNB),b(MNB,MNB)
      dimension pold(MNB,MNB),stemp(MNB,MNB)
      dimension nh_i(MNS),scr(MNB,MNB),iter_init(300)


      data nh_i/' S ',' P-',' P ',' D-',' D ',' F-',' F ',
     :                     ' G-',' G ','H-','H'/

c--------------------------------------------------------------------
c       initialize arrays and constants required for the computation.
c--------------------------------------------------------------------
      tcsqr=-two*c*c
      set = .true.
      iter=0.0
      en=zero
      eold=zero
      delta=zero

 
      do ia=1,MNB
      sa(ia)=zero
      do ib=1,MNB
      can(ib,ia)=zero
      drs(ib,ia)=zero
      p(ia,ib)=zero
      enddo
      enddo

      mbasis=2*nbasis
      write(STDOUT,199)maxit,npower
      write(STDOUT,200)    
      write(STDOUT,201)iter,delta
      write(*,201)iter,delta
 100  iter=iter+1
      do isym=1,nsym
      ndim=nbas(isym)
      ndim2=2*ndim
     

      do ia=1,ndim2
      ja=ia+2*nskip(isym)
      do ib=1,ndim2
      jb=ib+2*nskip(isym)
      stemp(ia,ib)=emt(ja,jb)
      enddo
      enddo
     
      if(jprint.ne.0)call matout(stemp,MNB,ndim2,ndim2,ndim2,4hS   )

      do i=1,ndim2
      do j=1,ndim2
      drs(i,j)=stemp(i,j)
      
      enddo
      enddo

      call eigen (drs,scr,ndim2,sa,ndim2,bnorm,MNB)
c-------------------------------------------------------------------
c       canonical transformation
c------------------------------------------------------------------
      

      do ia=1,ndim2
      do ib=1,ndim2
      if(dabs(sa(ib)).gt.TINY) then
      can(ia,ib)=drs(ia,ib)/dsqrt(dabs(sa(ib)))
      st(ib,ia)=can(ia,ib)
      else
      indt=1
      write(*,*)sa(ib)
      goto 705
      endif
      enddo
      enddo
      if(jprint.ne.0)call matout(st,MNB,ndim2,ndim2,ndim2,4hXT  )

      if(jprint.ne.0)
     :          call matout(df_single,MNB,mbasis,mbasis,mbasis,4hH   )

      if(iter.ne.1)call twoe(isym,nbasis,iter)


      do ic=1,ndim2
      jc=ic+2*nskip(isym)
      do id=1,ndim2
      jd=id+2*nskip(isym)
      if(iter.ne.1)then
      a(ic,id)=df_single(jc,jd)+df_two(jc,jd)
      else
      a(ic,id)=df_single(jc,jd)
      endif

      enddo
      enddo

      if(jprint.ne.0)call matout(a,MNB,ndim2,ndim2,ndim2,4hFOCK )

      if(jprint.ne.0)
     :call matout(df_two,MNB,mbasis,mbasis,mbasis,4hVHF )

      call mult(a,can,b,MNB,ndim2)
      call mult(st,b,df_mat,MNB,ndim2)
      if(jprint.ne.0)
     :      call matout(df_mat,MNB,ndim2,ndim2,ndim2,6hFRIME )
      call eigen (df_mat,scr,ndim2,eng,ndim2,bnorm,MNB)
      do i=1,ndim2
      do j=1,ndim2
      drs(i,j)=df_mat(i,j)
      enddo
      enddo

      if(jprint.ne.0)
     : call matout(drs,MNB,ndim2,ndim2,ndim2,4hC   )
      call mult(can,drs,eigv,MNB,ndim2)
  
      loc = 1
      last=ndim-1
      dowhile(loc.le.last)
      samin = eng(loc)
      locus = loc
      do  ia = loc,ndim
      if (eng(ia).lt.samin) then
      samin = eng(ia)
      locus = ia
      endif
      enddo
      if (locus .gt. loc) then
      swap = eng(loc)
      eng(loc) = eng(locus)
      eng(locus) = swap
      do  ib = 1,ndim2
      swap = eigv(ib,loc)
      eigv(ib,loc) = eigv(ib,locus)
      eigv(ib,locus) = swap
      enddo
      endif
      loc = loc + 1
      enddo
c-----------------------------------------------------------------
      do ia=1,ndim2
      do ib=1,ndim2
      b(ia,ib)=zero
      enddo
      b(ia,ia)=eng(ia)
      oe(ia,isym)=eng(ia)
      enddo

      if(jprint.ne.0)call matout(b,MNB,ndim2,ndim2,ndim2,4hE   )

      if(jprint.ne.0)call matout(eigv,MNB,ndim2,ndim2,ndim2,4hCP  )

      do ia=1,ndim2
      ja=ia+2*nskip(isym)
      do ib=1,ndim2
      jb=ib+2*nskip(isym)
      pold(ja,jb)=p(ja,jb)
      p(ja,jb)=0.0d0
      vec(ja,jb)=0.0d0
      nocc=nocorb(isym)
      if(nocc.ne.0)then
      do ibas=1,nocc
      twojp1=two*occ(isym)+one
      p(ja,jb)=p(ja,jb)+twojp1*eigv(ia,ibas)*eigv(ib,ibas)
      enddo
      endif
      vec(ja,jb)=eigv(ia,ib)
      enddo
      enddo

c   initialise 
       do ia =1,nbas(isym)
       ka=ia+nskip(isym)
       ja=ia+2*nskip(isym)
       do ic=1,n
       pf(ic,ka)=0.0d0
       qf(ic,ka)=0.0d0
       enddo
       ja=ia+2*nskip(isym)
       do ib =1,nbas(isym)
       jb=ib+2*nskip(isym)
       kb=ib+nskip(isym)
       do  ic=1,n
       pf(ic,ka)=pf(ic,ka)+vec(jb,ja)*gl(ic,kb)
       qf(ic,ka)=qf(ic,ka)+vec(jb+nbas(isym),ja)*gs(ic,kb)


       enddo
       enddo
       
       enddo

       enddo

      if(jprint.ne.0)call matout(p,MNB,mbasis,mbasis,mbasis,4hP   )

c     calculate delta

      delta=0.0d0
      do i=1,mbasis
      do j=1,mbasis
      delta=delta+(p(i,j)-pold(i,j))**2
      enddo
      enddo
c     indt=0
      delta=dsqrt(delta/4.0d0)
      write(STDOUT,201)iter,delta
      if(isk.eq.1)iter_init(isk)=iter
c     print*,'iter_init',iter_init
      if(isk.ne.1)iter_init(isk)=iter
      if(iter_init(1)+5.le.iter_init(isk))indt=1

      if(iter_init(1)+5.le.iter_init(isk))goto 705
      indt=0
        
      write(*,201)iter,delta
      
      if(delta.lt.crit)go to 160

      if(iter.lt.maxit)go to 100
      write(STDOUT,202)
      stop
  160 continue

      call integral(nbasis)
      eold=en
      en=zero
      do i=1,mbasis
      do j=1,mbasis
      en=en+0.5d0*p(i,j)*(two*df_single(i,j)+df_two(i,j))
      enddo
      enddo
      
      write(STDOUT,203)iter
      write(STDOUT,204)en
      write(STDOUT,205)
      do jsym=1,nsym
      do ibas=1,nbas(jsym)
      write(STDOUT,206)ibas,nh_i(jsym),oe(ibas,jsym),
     :                 oe(ibas+nbas(jsym),jsym)
       do i=1,itot
        if(i.le.ida(1))nsl=1
        if(i.le.ida(2).and.i.gt.ida(1))nsl=2
        if(i.le.ida(3).and.i.gt.ida(2))nsl=3
        if(i.le.ida(4).and.i.gt.ida(3))nsl=4
        if(i.le.ida(5).and.i.gt.ida(4))nsl=5
        if(i.le.ida(6).and.i.gt.ida(5))nsl=6
        if(i.le.ida(7).and.i.gt.ida(6))nsl=7
        if(i.le.ida(8).and.i.gt.ida(7))nsl=8
        if(i.le.ida(9).and.i.gt.ida(8))nsl=9
        if(i.le.ida(10).and.i.gt.ida(9))nsl=10
        if(i.le.ida(11).and.i.gt.ida(10))nsl=11
       if(ibas.eq.norb(i).and.jsym.eq.nsl) gauss1(i)= oe(ibas,jsym)
       enddo
      enddo
      enddo


 705     return
       
 199  format(/,1x,'CONVERGENCE DATA',/,
     1        ' MAXIMUM NO. OF ITERATION=',i6,/,
     2        ' CONVERGENCE CRITERION   = 1.0D-',i2,//)
 200  format(/,1x,'CYCLE      DENSITY CONV',/)
 201  format(i3,1x,(2x,d20.7))
 202  format(///,4x,33hSCF FAILS TO CONVERGES AT CYCLE , i4)
 203  format(///,4x,22hSCF CONVERGES AT CYCLE, i4)
 204  format(///,4x,20hELECTRONIC ENERGY=  ,d20.11)
 205  format(///,4x,32h ORBITAL ENERGIES (+ve and -ve) ,i4,///)
 206  format(/,4x,I2,A,2D20.11)
 207  format(/,4x,2D20.11)
      end            
************************************************************************
*                                                                      *
      double precision function rint (i,j,k)
*                                                                      *
*                                                                      *
*   the value of this rint is an approximation to:                     *
*                                                                      *
*              k                                                       *
*         i ( r  *  ( p (r)*p (r) + q (r)*q (r) ; 0 to infinity)       *
*                      i     j       i     j                           *
*                                                                      *
*   where   i ( g(r) ; range )  denotes  the  integral  of g(r) over   *
*   range.                                                             *
*                                                                      *
*   subroutines called: quad.                                          *
************************************************************************
*
      implicit real*8 (a-h,o-z)
c     include 'mpif.h'
      common/cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)
      dimension ta(MN)
*
      mtp = n
*
*   tabulate integrand as required for subroutine quad
*   value at first tabulation point is arbitrary
*
      ta(1) = zero
      do 1 l = 2,mtp
      ta(l)=(r(l)**k)*(pf(l,i)*pf(l,j)+qf(l,i)*qf(l,j))*rp(l)
    1 continue
*
*   perform integration
*
      call quad (ta,mtp,result)
      rint = result
*
      return
      end
c
      subroutine setgrd
      implicit real*8 (a-h,o-z)
c     include 'mpif.h'
      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n,lgrid
      common/quan1/kap(MNS),nak(MNOCC),np(MNOCC),nsy(MNS),nkk(MNS)


      do i = 2,MN
      r(i) = 0.0
      rp(i) = 0.0
      rpor(i) = 0.0
      enddo
      n = MN-10
      mtp=n
c     print*,'mtp',mtp
      h = 5.0d-02
      rnt = 2.0d-06
      r(1) = 0.0d0
      rp(1) = rnt
      eph = exp (h)
      ett = 1.d0
*   set up the arrays r, rp, rpor
      do i = 2,MN
      ett = eph*ett
      ettm1 = ett-1.0d0
      r(i) = rnt*ettm1
      rp(i) = rnt*ett
      rpor(i) = ett/ettm1
      enddo
      return
      end


      subroutine daclose
      implicit real*8 (a-h,o-z)
c     include 'mpif.h'

       common /double/df_two(MNB,MNB)
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c

      close (STDIN)
      close (STDOUT)
      close (IPUNCH)
      close (IWFN)
      close (15)
      close (19)
c     do i=1,MNB
c      do j=1,MNB
c     df_two(i,j)=0.0
c      enddo
c     enddo
      return
      end

