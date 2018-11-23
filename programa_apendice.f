      implicit real*8 (a-h,o-z)
      parameter(n=23,ns=7,md=n*n,neq=md)
      external fex, jex, base
      double complex rho(md),rhod(md),zwork(8*neq+2*neq*neq)
      double complex zhos(256)
      double precision rwork(10+neq)
      dimension nk(23,0:6), k(23), iwork(30+neq), pop(1:7,1:7)
      dimension qops(1:5,1:5), sops(1:5,1:5)
      double complex  rpar, wtru, err, az,ar,ai, t1,t2,t3,t4,ts    
      common /comp/ ar,ai,az
      common /const/ c, v,ggamma,omega0
      common /reserv/ gs,gd,th
      character gschar*6, gdchar*6, gammachar*5,vchar*8, cchar*5
      character omega0char*6, thchar*6
      open(80,file="entrada")
      read(80,*)c,v,gs,gd,tf,ggamma,th
      close(80)
      ar=(1d0,0d0)
      ai=(0d0,1d0)
      az=(0d0,0d0)
      pi=dacos(-1d0)
      rhoint5=0.d0
      nn=100
      call base(nk,k)  
      write(gammachar, "(i2.2,f4.3)") int(ggamma),ggamma-int(ggamma)
      write(gschar, "(i2.2,f4.3)") int(gs),gs-int(gs)
      write(gdchar, "(i2.2,f4.3)") int(gd),gd-int(gd)
      write(cchar, "(i2.2,f3.2)") int(c), c-int(c)
      write(vchar, "(i5.5,f3.2)") int(v), v-int(v)
      write(thchar, "(i2.2,f4.3)") int(th), th-int(th)
      open(51,file='sitio0_gs_'//gschar//'_gd_'//gdchar//'_v_'//vchar
     >     //'_th_'//thchar//'_gamma_'//gammachar//'.dat')
      open(52,file='sitio1_gs_'//gschar//'_gd_'//gdchar//'_v_'//vchar
     >     //'_th_'//thchar//'_gamma_'//gammachar//'.dat')
      open(53,file='sitio2_gs_'//gschar//'_gd_'//gdchar//'_v_'//vchar
     >     //'_th_'//thchar//'_gamma_'//gammachar//'.dat')
      open(54,file='sitio3_gs_'//gschar//'_gd_'//gdchar//'_v_'//vchar
     >     //'_th_'//thchar//'_gamma_'//gammachar//'.dat')
      open(55,file='sitio4_gs_'//gschar//'_gd_'//gdchar//'_v_'//vchar
     >     //'_th_'//thchar//'_gamma_'//gammachar//'.dat')
      open(56,file='sitio5_gs_'//gschar//'_gd_'//gdchar//'_v_'//vchar
     >     //'_th_'//thchar//'_gamma_'//gammachar//'.dat')
      open(57,file='sitio6_gs_'//gschar//'_gd_'//gdchar//'_v_'//vchar
     >     //'_th_'//thchar//'_gamma_'//gammachar//'.dat')
      open(58,file='media_gs_'//gschar//'_gd_'//gdchar//'_v_'//vchar
     >     //'_th_'//thchar//'_gamma_'//gammachar//'.dat')
      open(59,file='variancia_gs_'//gschar//'_gd_'//gdchar//'_v_'//
     >     vchar //'_th_'//thchar//'_gamma_'//gammachar//'.dat')
      open(60,file='traco2.dat')
      open(61, file='traco1.dat')
      open(62, file='soma_termos_cruzados.dat')
      open(63, file='traco_sistema.dat')
      kin=1
      t=0.d0
      do i=1,md
         rho(i)=az
      enddo
      do iis=1,256 
         zhos(iis)=az
      enddo
      do ik=1,n
         do ikl=1,n
            ii=(ik-1)*n+ikl
            delta=1d0
            delta1=1d0
            do m=0,ns-1
               adif=dabs(dfloat(nk(ik,m)-nk(kin,m)))
               if(adif.lt.0.0001d0) then
                  adif=1d0
               else
                  adif=0d0
               endif
               delta=adif*delta
               adif1=dabs(dfloat(nk(ikl,m)-nk(kin,m)))
               if(adif1.lt.0.0001d0) then
                  adif1=1d0
               else
                  adif1=0d0
               endif
               delta1=adif1*delta1
            enddo
            rho(ii)=rho(ii)+ar*delta*delta1
         enddo
      enddo
      soma=0.d0
      do i=1,n
         soma=soma+dreal(rho((i-1)*n+i))
      enddo
      write(8,*) 'traco inicial = ', soma
      dtout = 0.01d0
      tout = t+dtout
      itol = 1
      rtol = 1.d-11
      atol = 1.d-10
      itask = 1
      istate = 1
      iopt = 0
      lzw = 8*neq+2*neq*neq
      lrw = 20+neq
      liw = 30+neq
      mf = 21                   ! da rotina.
      rpar = dcmplx(0.0d0,1.0d0)
      aemax = 0.0d0
      loopmax=int((tf-t)/dtout)
      do 40 iout = 1,loopmax/nn
         do 60 j=1,nn
            call zvode(fex,neq,rho,t,tout,itol,rtol,atol,itask,istate,io
     >           pt,
     1           zwork,lzw,rwork,lrw,iwork,liw,jex,mf,rpar,ipar)
            wtru = 1.0d0/dcmplx(cos(t) + 1.1d0, sin(t))
            err = rho(1) - wtru
            aberr = abs(dreal(err)) + abs(dimag(err))
            aemax = max(aemax,aberr)
            do i=0,ns-1
               do  jj=0,ns-1
                  pop(i+1,jj+1)=0l.d0
                  if (i.eq.jj) then
                     do ik=1,n
                        iii=(ik-1)*n+ik
                        pop(i+1,jj+1)=pop(i+1,jj+1)+nk(ik,i)
     >                       *dreal(rho(iii))
                     end do
                  end if   
                  if (i.ne.jj) then
                     do ik=1,n
                        do ikl=1,n
                           if (ik.ne.ikl) then  
                              iii=(ik-1)*n+ikl
                              pop(i+1,jj+1)=pop(i+1,jj+1)+nk(ik,i)*
     >                             nk(ikl,jj)*abs(rho(iii))
                                end if 
                        end do
                     enddo
                  endif
               enddo
            enddo
            zhos(1)=rho(1)+rho(385)+rho(529)
            zhos(18)=rho(25)+rho(409)
            zhos(19)=rho(26)+rho(410)
            zhos(21)=rho(28)+rho(411)      
            zhos(24)=rho(31)+rho(412)
            zhos(28)=rho(35)+rho(413)
            zhos(34)=rho(48)+rho(432)
            zhos(35)=rho(49)+rho(433)
            zhos(37)=rho(51)+rho(434)
            zhos(40)=rho(54)+rho(435)
            zhos(44)=rho(58)+rho(436)
            zhos(52)=rho(73)
            zhos(54)=rho(75)
            zhos(55)=rho(76)
            zhos(57)=rho(78)
            zhos(58)=rho(79)
            zhos(59)=rho(80)
            zhos(61)=rho(82)
            zhos(62)=rho(83)
            zhos(63)=rho(84)
            zhos(64)=rho(85)
            zhos(66)=rho(94)+rho(455)
            zhos(67)=rho(95)+rho(456)
            zhos(69)=rho(97)+rho(457)
            zhos(72)=rho(100)+rho(458)
            zhos(76)=rho(104)+rho(459)
            zhos(84)=rho(119)
            zhos(86)=rho(121)
            zhos(87)=rho(122)
            zhos(89)=rho(124)
            zhos(90)=rho(125)
            zhos(91)=rho(126)
            zhos(93)=rho(128)
            zhos(94)=rho(129)
            zhos(95)=rho(130)
            zhos(96)=rho(131)
            zhos(100)=rho(142)
            zhos(102)=rho(144)
            zhos(103)=rho(145)
            zhos(105)=rho(147)
            zhos(106)=rho(148)
            zhos(107)=rho(149)
            zhos(109)=rho(151)
            zhos(110)=rho(152)
            zhos(111)=rho(153)
            zhos(112)=rho(154)
            zhos(114)=rho(163)+rho(478)
            zhos(115)=rho(164)+rho(479)
            zhos(117)=rho(166)+rho(480)
            zhos(120)=rho(169)+rho(481)
            zhos(124)=rho(173)+rho(482)
            zhos(132)=rho(188)
            zhos(134)=rho(190)
            zhos(135)=rho(191)
            zhos(137)=rho(193)
            zhos(138)=rho(194)
            zhos(139)=rho(195)
            zhos(141)=rho(197)
            zhos(142)=rho(198)
            zhos(143)=rho(199)
            zhos(144)=rho(200)
            zhos(148)=rho(211)
            zhos(150)=rho(213)
            zhos(151)=rho(214)
            zhos(153)=rho(216)
            zhos(154)=rho(217)
            zhos(155)=rho(218)
            zhos(157)=rho(220)
            zhos(158)=rho(221)
            zhos(159)=rho(222)
            zhos(160)=rho(223)
            zhos(164)=rho(234)
            zhos(166)=rho(236)
            zhos(167)=rho(237)
            zhos(169)=rho(239)
            zhos(170)=rho(240)
            zhos(171)=rho(241)
            zhos(173)=rho(243)
            zhos(174)=rho(244)
            zhos(175)=rho(245)
            zhos(176)=rho(246)
            zhos(178)=rho(255)+rho(501)
            zhos(179)=rho(256)+rho(502)
            zhos(181)=rho(258)+rho(503)
            zhos(184)=rho(261)+rho(504)
            zhos(188)=rho(265)+rho(505)
            zhos(196)=rho(280)
            zhos(198)=rho(282)
            zhos(199)=rho(283)
            zhos(201)=rho(285)
            zhos(202)=rho(286)
            zhos(203)=rho(287)
            zhos(205)=rho(289)
            zhos(206)=rho(290)
            zhos(207)=rho(291)
            zhos(208)=rho(292)
            zhos(212)=rho(303)
            zhos(214)=rho(305)
            zhos(215)=rho(306)
            zhos(217)=rho(308)
            zhos(218)=rho(309)
            zhos(219)=rho(310)
            zhos(221)=rho(312)
            zhos(222)=rho(313)
            zhos(223)=rho(314)
            zhos(224)=rho(315)
            zhos(228)=rho(326)
            zhos(230)=rho(328)
            zhos(231)=rho(329)
            zhos(233)=rho(331)
            zhos(234)=rho(332)
            zhos(235)=rho(333)
            zhos(237)=rho(335)
            zhos(238)=rho(336)
            zhos(239)=rho(337)
            zhos(240)=rho(338)
            zhos(244)=rho(349)
            zhos(246)=rho(351)
            zhos(247)=rho(352)
            zhos(249)=rho(354)
            zhos(250)=rho(355)
            zhos(251)=rho(356)
            zhos(253)=rho(358)
            zhos(254)=rho(359)
            zhos(255)=rho(360)
            zhos(256)=rho(361)
            pop_t=0d0
            do i=1,7
               pop_t=pop_t+pop(i,i)
            enddo
            rhoint5=rhoint5+gd*pop(6,6)*dtout
            vmedia=0.d0
            vmedia2=0.d0
            do i=1,7
               vmedia=vmedia+(i)*pop(i,i)/2.d0
               vmedia2=vmedia2+(i**2)*pop(i,i)/2.d0
            enddo   
            do is=1,5
               do  jjs=1,5
                  qops(is,jjs)=0.d0
                  sops(is,jjs)=0.d0
                  if (is.eq.jjs) then
                     do iks=1,16
                        iiis=(iks-1)*16+iks
                        qops(is,jjs)=qops(is,jjs)+nk(iks,is)
     >                       *dreal(zhos(iiis))
                     end do
                  end if   
                  if (is.ne.jjs) then
                     do iks=1,16
                        do ikls=1,16
                           if (iks.ne.ikls) then  
                              iiis=(iks-1)*16+ikls
                              qops(is,jjs)=qops(is,jjs)+nk(iks,is)*
     >                             nk(ikls,jjs)*abs(zhos(iiis))
                           end if 
                        end do
                     enddo
                  endif
               enddo
            enddo    
 60         tout = tout + dtout      
            soma=0.d0
            do i=1,n
               soma=soma+dreal(rho((i-1)*n+i))
            enddo
            pops_t=0d0
            do i=1,5
               pops_t=pops_t+qops(i,i)
            enddo  
            do i=1,5
               do ji=1,5
                  if (i.eq.ji) then
                     sops(i,ji)=0.d0
                  else
                     sops(i,ji)=qops(i,ji)
                  endif         
               enddo
               write(82, '(5(f10.8,1x))') (sops(i,j), j=1,5)
            enddo
            write(82, *) t
            ssops=0.d0
            do in=1,5
               do jn=1,5
                  ssops=ssops+sops(in,jn)
               enddo
            enddo
            write(51,*) t, pop(1,1)
            write(52,*) t, pop(2,2)
            write(53,*) t, pop(3,3)
            write(54,*) t, pop(4,4)
            write(55,*) t, pop(5,5)
            write(56,*) t, pop(6,6)
            write(57,*) t, pop(7,7)
            write(58,*) t, vmedia
            write(59,*) t, (vmedia2)-(vmedia)**2
            write(60,*) t, pop_t
            write(61,*) t, soma
            write(62,*) t, ssops
            write(63,*) t, pops_t
            do i=1,n
               write(72,*) t, rho((i-1)*n+i)
            enddo
            do i=1,16
               write(92,*) t, zhos((i-1)*16+i)
            enddo
            call flush()        
 20         format(f5.1,2x,2f9.5,3x,2f9.5)
 40      continue
         stop
         end
      subroutine fex (neq, t, rho, rhod, rpar, ipar)
      implicit real*8 (a-h,o-z)
      double complex rho(neq), rhod(neq), rpar(1)
      double complex ar,ai,az, t1,t2,t3,t4,ts   
      external base
      dimension nk(23,0:6), k(23)
      integer ipar(1)
      common /comp/ ar,ai,az
      common /const/ c, v,ggamma, omega0
      common /reserv/ gs,gd,th
      n=23
      ns=7
      call base(nk,k)
      do ik=1,n
         do ikl=1,n
            ii=(ik-1)*n+ikl
            rhod(ii)=az
         enddo
      enddo
      do ik=1,n
         do ikl=1,n
            ii=(ik-1)*n+ikl
c.....abre termo do campo
            do jj=1,6
               rhod(ii)=rhod(ii)+ai*(omega0)*dfloat(jj)*(nk(ik,jj)-
     >              nk(ikl,jj))*rho(ii)
            enddo
c.....fecha termo do campo    
c.....abre termo de tunelamento
            do j=1,4            !nenhum tunelamento para fonte ou sumidouro
               it1=k(ik)+2*3ll**j
               do m=1,n
                  km=k(m)
                  x=dabs(dfloat(it1)-dfloat(km))
                  if(x.le.0.001d0) then
                     im=(m-1)*n+ikl
                     t1=rho(im)
                     go to 1
                  else
                     t1=az
                  endif
               enddo
 1             continue
               it2=k(ik)-2*3**j
               do m=1,n
                  km=k(m)
                  x=dabs(dfloat(it2)-dfloat(km))
                  if(x.le.0.001d0) then
                     im=(m-1)*n+ikl
                     t2=rho(im)
                     go to 2
                  else
                     t2=az
                  endif
               enddo
 2             continue
               it3=k(ikl)-2*3**j
               do m=1,n
                  km=k(m)
                  x=dabs(dfloat(it3)-dfloat(km))
                  if(x.le.0.001d0) then
                     im=(ik-1)*n+m
                     t3=rho(im)
                     go to 3
                  else
                     t3=az
                  endif
               enddo
 3             continue
               it4=k(ikl)+2*3**j
               do m=1,n
                  km=k(m)
                  x=dabs(dfloat(it4)-dfloat(km))
                  if(x.le.0.001d0) then
                     im=(ik-1)*n+m
                     t4=rho(im)
                     go to 4
                  else
                     t4=az
                  endif
               enddo
 4             continue
               rhod(ii)=rhod(ii)+ai*c*(nk(ik,j)*(1-nk(ik,j+1))*t1
     >              +nk(ik,j+1)*(1-nk(ik,j))*t2
     >              -(1-nk(ikl,j))*nk(ikl,j+1)*t3
     >              -nk(ikl,j)*(1-nk(ikl,j+1))*t4)
            enddo
c.....fecha termo de tunelamento
c.....abre termo de coulomb 
            do j=1,4
               do jl=j+1,5
                  a=abs(j-jl)
                  
                  rhod(ii)=rhod(ii)+ai*(v/a)*(nk(ik,j)*nk(ik,jl)-
     >                 nk(ikl,j)*nk(ikl,jl))*rho(ii)
               enddo
            enddo
c......fecha termo de coulomb
c-----abre fonte
            its1=k(ik)-2
            its2=k(ikl)-2
            if(its1.ge.2.and.its1.le.1458) then
               if(its2.ge.2.and.its2.le.1458) then
                  do m=1,n
                     km=k(m)
                     
                     if(its1.eq.km) then
                        im1=m
                        aa=dfloat(1+nk(ik,0))
                        aux1=dsqrt(aa)*nk(ik,1)
                        go to 50
                     else
                        im1=0
                     endif
                  enddo
 50               continue
                  do m=1,n
                     km=k(m)
                     
                     if(its2.eq.km) then
                        im2=m
                        go to 51
                     else
                        im2=0
                     endif
                  enddo
 51               continue
                  iprod=im1*im2
                  if(iprod.gt.0) then
                     im=(im1-1)*n+im2
                     ts=rho(im)
                     aa=dfloat(1+nk(ikl,0))
                     rhod(ii)=rhod(ii)+2.d0*gs*aux1*nk(ikl,1)*dsqrt(aa)*ts
                     go to 52
                  endif
               endif
            else
               ts=az
            endif
 52   continue
      rhod(ii)=rhod(ii)-gs*dfloat(nk(ik,0)*(1-nk(ik,1))+nk(ikl,0)*
     >     (1-nk(ikl,1)))*rho(ii)
c-----fecha fonte
c-----abre sumidouro
      itd1=k(ik)-486
      itd2=k(ikl)-486
      if(itd1.ge.2.and.itd1.le.1458) then
         if(itd2.ge.2.and.itd2.le.1458) then
            do m=1,n
               km=k(m)
               if(itd1.eq.km) then
                  im1=m
                  aa=dfloat(nk(ik,6))
                  aux1=dsqrt(aa)*(1-nk(ik,5))
                  go to 62
               else
                  im1=0
               endif
            enddo
 62         continue
            do m=1,n
               km=k(m)
               if(itd2.eq.km) then
                  im2=m
                  go to 63
               else
                  im2=0
               endif
            enddo
 63         continue
            iprod=im1*im2
            if(iprod.gt.0) then
               im=(im1-1)*n+im2
               ts=rho(im)
               aa=dfloat(nk(ikl,6))
               rhod(ii)=rhod(ii)+2.d0*gd*aux1*(1-nk(ikl,5))*dsqrt(aa)*ts
               go to 222
            endif
         endif
      endif
 222  continue
      rhod(ii)=rhod(ii)-gd*dfloat(nk(ik,5)*(1+nk(ik,6))+
     >     nk(ikl,5)*(1+nk(ikl,6)))*rho(ii)   
c-----fecha sumidouro
c-----abre dephasing
      do jj=1,5
         rhod(ii)=rhod(ii)+ggamma*(2*nk(ik,jj)*nk(ikl,jj)-nk(ik,jj)
     >        -nk(ikl,jj))*rho(ii)
      enddo
c-----fecha dephasing 
c-----abre thermal
      do j=1,4 
         itth1=k(ik)-2*(3**j)
         itth2=k(ikl)-2*(3**j)
         do m=1,n
            km=k(m)
            if(itth1.eq.km) then
               imth1=m
               go to 11
            else
               imth1=0
            endif
         enddo
 11      continue
         do m=1,n
            km=k(m)
            if(itth2.eq.km) then
               imth2=m
               go to 22
            else
               imth2=0
            endif
         enddo
 22      continue
         iprodth=imth1*imth2
         if(iprodth.gt.0)then
            imth=(imth1-1)*n+imth2
            tt1=rho(imth)
         end if
         itth3=k(ik)+2*(3**j)
         itth4=k(ikl)+2*(3**j)
         do m=1,n
            km=k(m)
            if(itth3.eq.km) then
               imth13=m
               go to 33
            else
               imth13=0
            endif
         enddo
 33      continue
         do m=1,n
             km=k(m)
             if(itth4.eq.km) then
                imth14=m
                go to 44
             else
                imth14=0
             endif
          enddo
 44       continue
          iiprodth=imth13*imth14
          if(iiprodth.gt.0)then
             imth=(imth13-1)*n+imth14
             tt2=rho(imth)
          end if
          rhod(ii)=rhod(ii)+th*(2*nk(ik,j+1)*(1-nk(ik,j))
     >         *nk(ikl,j+1)*(1-nk(ikl,j))*tt1
     >         -((1-nk(ik,j+1))*nk(ik,j)
     >         +nk(ikl,j)*(1-nk(ikl,j+1)))*rho(ii))
     >         +th*(2*nk(ik,j)*(1-nk(ik,j+1))
     >         *nk(ikl,j)*(1-nk(ikl,j+1))*tt2
     >         -((1-nk(ik,j))*nk(ik,j+1)
     >         +nk(ikl,j+1)*(1-nk(ikl,j)))*rho(ii))
       enddo
c---- fecha thermal 
      enddo
      enddo
      return
      end
      subroutine jex (neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
      implicit real*8 (a-h,o-z)
      double complex  y(neq), pd(nrpd,neq), rpar(1)
      dimension nk(23,0:6), k(23) 
      integer ipar(1)          
      double complex ar,ai,az, t1,t2,t3,t4,ts
      external base
      common /comp/ ar,ai,az
      common /const/ c, v,ggamma, omega0
      common /reserv/ gs,gd, th
      n=23
      ns=7
      call base(nk,k)
      do i=1,neq
         do j=1,neq
            pd(i,j)=az
         enddo
      enddo
c.....abre termo de tunelamento
      do ik=1,n
         do ikl=1,n
            ii=(ik-1)*n+ikl
c.....abre termo do campo
            do jj=1,6
               pd(ii,ii)=pd(ii,ii)+ai*omega0*dfloat(jj)*(nk(ik,jj)-
     >              nk(ikl,jj))
            enddo
c.....fecha termo do campo  
        do j=1,4 !sem tunelamento para fonte ou sumidouro
         it1=k(ik)+2*3**j
          do m=1,n
           km=k(m)
           x=dabs(dfloat(it1)-dfloat(km))
           if(x.le.0.001d0) then
            im=(m-1)*n+ikl
            im1=im
            t1=ar
            go to 1
           else
            t1=az
           endif
          enddo
 1       continue
         it2=k(ik)-2*3**j
          do m=1,n
           km=k(m)
           x=dabs(dfloat(it2)-dfloat(km))
           if(x.le.0.001d0) then
             im=(m-1)*n+ikl
             im2=im
             t2=ar
             go to 2
            else
             t2=az
           endif
          enddo
 2       continue
         it3=k(ikl)-2*3**j
          do m=1,n
           km=k(m)
           x=dabs(dfloat(it3)-dfloat(km))
           if(x.le.0.001d0) then
             im=(ik-1)*n+m
             im3=im
             t3=ar
             go to 3
            else
             t3=az
           endif
          enddo
 3       continue
         it4=k(ikl)+2*3**j
          do m=1,n
           km=k(m)
           x=dabs(dfloat(it4)-dfloat(km))
           if(x.le.0.001d0) then
             im=(ik-1)*n+m
             im4=im
             t4=ar
             go to 4
            else
             t4=az
           endif
          enddo
 4       continue
         pd(ii,im1)=pd(ii,im1)+ai*c*nk(ik,j)*(1-nk(ik,j+1))*t1
         pd(ii,im2)=pd(ii,im2)+ai*c*nk(ik,j+1)*(1-nk(ik,j))*t2
         pd(ii,im3)=pd(ii,im3)-ai*c*(1-nk(ikl,j))*nk(ikl,j+1)*t3
         pd(ii,im4)=pd(ii,im4)-ai*c*nk(ikl,j)*(1-nk(ikl,j+1))*t4
        enddo
c.....fecha termo de tunelamento
c.....abre termo de coulomb 
         do j=1,4
          do jl=j+1,5
          a=abs(j-jl)
           pd(ii,ii)=pd(ii,ii)+ai*(v/a)*(nk(ik,j)*nk(ik,jl)-
     >      nk(ikl,j)*nk(ikl,jl))
         enddo
        enddo
c.....fecha termo de coulomb
c---- abre termo de fonte
        its1=k(ik)-2
        its2=k(ikl)-2
        if(its1.ge.2.and.its1.le.1458) then
           if(its2.ge.2.and.its2.le.1458) then
              do m=1,n
                 km=k(m)
                 if(its1.eq.km) then
                    im1=m
                    aa=dfloat(1+nk(ik,0))
                    aux1=dsqrt(aa)*nk(ik,1)
                    go to 81
                 else
                    im1=0 
                 endif
              enddo
 81           continue
              do m=1,n
                 km=k(m)
                 if(its2.eq.km) then
                    im2=m
                    go to 82
                 else
                    im2=0
                 endif
              enddo
 82           continue
              iprod=im1*im2
              if(iprod.gt.0) then
                 im=(im1-1)*n+im2
                 ts=ar
                 aa=dfloat(1+nk(ikl,0))
                 pd(ii,im)=pd(ii,im)+2.d0*gs*aux1*nk(ikl,1)*dsqrt(aa)*ts
                 go to 881
              endif
           endif
        else
           ts=az
        endif
 881    continue
        pd(ii,ii)=pd(ii,ii)-gs*dfloat(nk(ik,0)*(1-nk(ik,1))+nk(ikl,0)*
     >       (1-nk(ikl,1)))*ar
c-----fecha fonte
c-----abre termo de sumidouro 
        itd1=k(ik)-486
        itd2=k(ikl)-486
        if(itd1.ge.2.and.itd1.le.1458) then
           if(itd2.ge.2.and.itd2.le.1458) then
              do m=1,n
                 km=k(m)
                 if(itd1.eq.km) then
                    im1=m
                    aa=dfloat(nk(ik,6))
                    aux1=dsqrt(aa)*(1-nk(ik,5))
                    go to 83
                 else
                    im1=0 
                 endif
              enddo
 83           continue
              do m=1,n
                 km=k(m)
           if(itd2.eq.km) then
              im2=m
              go to 84
           else
              im2=0
           endif
        enddo
 84     continue
        iprod=im1*im2
        if(iprod.gt.0) then
           im=(im1-1)*n+im2
           ts=ar
           aa=dfloat(nk(ikl,6))
           pd(ii,im)=pd(ii,im)+2d0*gd*aux1*(1-nk(ikl,5))*
     >          dsqrt(aa)*ts
           go to 99
        endif
         endif
        endif
 99     continue
        pd(ii,ii)=pd(ii,ii)-gd*dfloat(nk(ik,5)*(1+nk(ik,6))+
     >   nk(ikl,5)*(1+nk(ikl,6)))*ar
c-----fecha sumidouro
c--------abre dephasing
        do jj=1,5
           pd(ii,ii)=pd(ii,ii)+ggamma*(2*nk(ik,jj)*nk(ikl,jj)-nk(ik,jj)
     v          -nk(ikl,jj))*ar
        enddo
c-------fecha dephasing        
c-------abre thermal
        do j=1,4 
           itth1=k(ik)-2*3**j
           itth2=k(ikl)-2*3**j
           do m=1,n
              km=k(m)
              if(itth1.eq.km) then
                 imth1=m
                 go to 11
              else
                 imth1=0
              endif
           enddo
 11        continue
           do m=1,n
              km=k(m)
              if(itth2.eq.km) then
                 imth2=m
                 go to 22
              else
                 imth2=0
              endif
           enddo
 22        continue
           iprodth=imth1*imth2
           if(iprodth.gt.0)then
              ts=ar
              imth=(imth1-1)*n+imth2
              pd(ii,imth)=pd(ii,imth)+2*th*nk(ik,j+1)*(1-nk(ik,j))
     >                               *nk(ikl,j+1)*(1-nk(ikl,j))*ts
              go to 101
           end if
 101       pd(ii,ii)=pd(ii,ii)-th*((1-nk(ik,j+1))*nk(ik,j)
     >          +nk(ikl,j)*(1-nk(ikl,j+1)))*ar
           itth3=k(ik)+2*3**j
           itth4=k(ikl)+2*3**j
           do m=1,n
              km=k(m)
              if(itth3.eq.km) then
                 imth13=m
                 go to 33
              else
                 imth13=0
              endif
           enddo
 33        continue
           do m=1,n
              km=k(m)
              if(itth4.eq.km) then
                 imth14=m
                 go to 44
              else
                 imth14=0
              endif
           enddo
 44        continue
           iiprodth=imth13*imth14
           if(iiprodth.gt.0)then
              ts=ar
              imth=(imth13-1)*n+imth14
              pd(ii,imth)=pd(ii,imth)+2*th*nk(ik,j)*(1-nk(ik,j+1))
     >                               *nk(ikl,j)*(1-nk(ikl,j+1))*ts
              go to 104
           end if
 104       continue
           pd(ii,ii)=pd(ii,ii)-th*((1-nk(ik,j))*nk(ik,j+1)
     >          +nk(ikl,j+1)*(1-nk(ikl,j)))*ar
        enddo
c-------fecha thermal        
       enddo
      enddo
      return
      end
      subroutine base(nk,k)
      implicit real*8 (a-h,o-z)
      dimension nk(23,0:6), k(23)
      do j=1,23
       do i=0,6
        nk(j,i)=0
       enddo
      enddo
      nk(1,0)=2
      nk(2,0)=1
      nk(2,1)=1
      nk(3,0)=1
      nk(3,2)=1
      nk(4,1)=1
      nk(4,2)=1
      nk(5,0)=1
      nk(5,3)=1
      nk(6,1)=1
      nk(6,3)=1
      nk(7,2)=1
      nk(7,3)=1
      nk(8,0)=1
      nk(8,4)=1
      nk(9,1)=1
      nk(9,4)=1
      nk(10,2)=1
      nk(10,4)=1
      nk(11,3)=1
      nk(11,4)=1
      nk(12,0)=1
      nk(12,5)=1
      nk(13,1)=1
      nk(13,5)=1
      nk(14,2)=1
      nk(14,5)=1
      nk(15,3)=1
      nk(15,5)=1
      nk(16,4)=1
      nk(16,5)=1
      nk(17,0)=1
      nk(17,6)=1
      nk(18,1)=1
      nk(18,6)=1
      nk(19,2)=1
      nk(19,6)=1
      nk(20,3)=1
      nk(20,6)=1
      nk(21,4)=1
      nk(21,6)=1
      nk(22,5)=1
      nk(22,6)=1
      nk(23,6)=2
      do ik=1,23
       k(ik)=0
       do j=0,6
        k(ik)=nk(ik,j)*3**j+k(ik)
       enddo
      enddo
      return
      end

