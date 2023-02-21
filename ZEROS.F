C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (10/07/2009)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c                     initial values specifying                             #
c                                                                           #
c############################################################################

      Subroutine ZEROES
	
      Include './Include/OCERM_INF'
      
      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,XO,XL,XR,FYI)
!$OMP DO
      Do I = 1, IJM
         Do K = 1, KB
            T(I,K) = 0.0
            SED(I,K) = 0.0
	      SAL(I,K) = 0.0
	      WQ(I,K) = 0.0 

C            P1(I,K)=0.0
C            Z1(I,K)=0.0
C            N1(I,k)=0.0

C            RHO(I,K) = 0.0
            RMEAN(I,K) = 1000.0
C            TMEAN(I,J,K) = 0.0
C            SMEAN(I,J,K) = 0.0
C	      RI(I,J,K) = 0.0

C            PMEAN(I,J,K) =0.0
C            ZMEAN(I,J,K) =0.0
C            NMEAN(I,J,K) =0.0

C            A(I,J,K) = 0.0
C            C(I,J,K) = 0.0
C            VH(I,J,K) = 0.0
C            VHP(I,J,K) = 0.0
C            PROD(I,J,K) = 0.0
C            DTEF(I,J,K) = 0.0
            U(I,K) = 0.0
            V(I,K) = 0.0
	      W(I,K) = 0.0
	      QZ(I,K) = 0.0
            UR(I,K) = 0.0
            VR(I,K) = 0.0
	      WR(I,K) = 0.0
            UF(I,K) = 0.0
            VF(I,K) = 0.0
!----------#Turbulence tolerance#, WangJian, 2020-3-7 12:08:09-----------------£¡	      
            TKE(I,K) = 1.0E-20
	      TDISS(I,K) = 1000.0           
            VIS(I,K) = 1.0E-10
C	      Do N = 1, 3
C	         Do M = 1, 3
C	            STRESS(M,N,I,K) = 0.0
C	         Enddo
C	      Enddo
C            UBETA(I,J,K) = 0.0
C            VBETA(I,J,K) = 0.0
            
C            DRHOX(I,J,K) = 0.0
C            DRHOY(I,J,K) = 0.0
C            Q2(I,J,K) = 0.0
C            Q2L(I,J,K) = 0.0
C            L(I,J,K) = 0.0
            KM(I,K) = 0.0
            KH(I,K) = 0.0
            KQ(I,K) = 0.0

            DELP(I,K,1) = 0.0
            DELP(I,K,2) = 0.0
            DELP(I,K,3) = 0.0           
            
	      GRADP(I,K,1) = 0.0
	      GRADP(I,K,2) = 0.0 
            GRADP(I,K,3) = 0.0             
            
C	      XMFLUX(I,J,K)=0.0
C	      YMFLUX(I,J,K)=0.0

C	      RADSX(I,J,K)=0.0
C	      RADSY(I,J,K)=0.0
        Enddo
      Enddo
!$OMP END DO NOWAIT
!$OMP DO
	Do I = 1, IJE
	   Do K = 1, KB
	      UN(I,K) = 0.0
	      VN(I,K) = 0.0
	   Enddo
	Enddo
!$OMP END DO NOWAIT
!$OMP DO
	Do I = 1, IJP
	   Do K = 1, KB
	      UV(I,K) = 0.0
	      VV(I,K) = 0.0
	      WV(I,K) = 0.0
	   Enddo
	Enddo
!$OMP END DO NOWAIT
C---------------------------------------------------------------------------------C
C         make solitary wave                                                      C
C---------------------------------------------------------------------------------C
!$OMP DO
	Do I = 1, IJM
         ELF(I) = 0.0
         EL(I) = 0.0
         X0 = 2.0   ! ²¨·åÎ»ÖÃ
         XL = 10.0    ! ²¨×ó¶Ë
         XR = 4.0    ! ²¨ÓÒ¶Ë
         HDEPTH = 0.75   ! Ë®Éî
         HWAVE = 0.3     ! ¹ÂÁ¢²¨²¨¸ß
         IF(CXY(I,1) .GE. XL .AND. CXY(I,1) .LE. XR) THEN
            FYI = SQRT(3.*HWAVE/(4.*HDEPTH**3.)) * (CXY(I,1)-X0)
            SS = 1. / COSH(FYI)
            TT = TANH(FYI)
     
            AMP_C = 0.13158   ! Ïà¶Ô²¨¸ß
            ELF(I) = HDEPTH*(AMP_C * SS**2. - 
     &               AMP_C**2. * 0.75 * SS**2. * TT**2. + 
     &           AMP_C**3.*(0.625*SS**2.*TT**2.-101./80.*SS**4.*TT**2.))
     
            EL(I) = ELF(I)
            Do K = 1, KBM
                ZK = 1. + ZZ(K)
               U(I,K) = SQRT(9.8*HDEPTH)*(HC(I) + EL(I)) *
     &                  (AMP_C * SS ** 2.0 + 
     &                   AMP_C**2. * (-0.25*SS**2.+SS**2.*TT**2.+
     &               ZK**2.*(1.5*SS**2.-2.25*SS**2.*TT**2.))+
     &  AMP_C**3. * (21.0/40.0*SS**2.-SS**2.*TT**2.-
     &  6./5.*SS**4.*TT**2.+
     &    ZK**2.*(-1.5*SS**2.+3.75*SS**2.*TT**2.+7.5*SS**4.*TT**2.)+
     &     ZK**4.*(0.375*SS**2.+45./16.*SS**4.*TT**2.)))

               QZ(I,K) = SQRT(9.8) * (HC(I) + EL(I)) *
     &                  SQRT(3.*AMP_C)*(1.+Z(K))*SS **2. * TT * 
     &              (AMP_C+0.*AMP_C**2.*(-0.375-2.*SS**2.+(1.+Z(K))**2.*
     &              (-0.5+1.5*SS**2.))+0.*
     &              AMP_C**3.*(-49./640.-17./20.*SS**2.-18./5.*SS**4.+
     &           ZK**2.*(-13./16.+25./16.*SS**2.+15./2.*SS**4.)+
     &           ZK**4.*(-3./40.+9./8.*SS**2.-27./16.*SS**4.)))
           ENDDO
            
         ENDIF
C         ELF(I) = 0.02*COS(31.4*CXY(I,1))
C        TPS(I) = 0.0
C	   FLUX(I)=0.0
         WTSURF(I,1) = 0.0
         WTSURF(I,2) = 0.0
C        SWRAD(I,J) = 0.0
C        WSSURF(I,J) = 0.0
         WUSURF(I) = 0.0
         WVSURF(I) = 0.0

	   DHT(I) = 0.0
	   IWB(I) = 0
C        WPSURF(I,J) = 0.0
C        WZSURF(I,J) = 0.0
C        WNSURF(I,J) = 0.0

C        WUBOT(I,J) = 0.0
C        WVBOT(I,J) = 0.0

C	   WA(I,J)=0.00
C	   WH(I,J)=0.00
C	   WNX(I,J)=0.0
C	   WNY(I,J)=0.0
C	   WFE(I,J)=0.0

         Z01(I)=0.0
C	   H0(I,J)=0.0
C	   HT(I,J)=0.0
      Enddo
!$OMP END DO
!$OMP END PARALLEL

	Do I = 1, IJM
	   DC(I) = HC(I) + EL(I)
	Enddo
	Do I = 1, IJE
	   If(INDEX_EDGE(I,1,1) .EQ. -999) Then
       	  DS(I) = HS(I) + EL(INDEX_EDGE(I,1,2))
	      ELFM(I) = EL(INDEX_EDGE(I,1,2))
	   Endif
	   If(INDEX_EDGE(I,1,2) .EQ. -999) Then
      	  DS(I) = HS(I) + EL(INDEX_EDGE(I,1,1))
	      ELFM(I) = EL(INDEX_EDGE(I,1,1))
	   Endif
	   If(INDEX_EDGE(I,1,1).NE.-999 .AND. INDEX_EDGE(I,1,2).NE.-999) 
     &   Then
            DS(I) = HS(I) + 0.5 * 
     &              (EL(INDEX_EDGE(I,1,1)) + EL(INDEX_EDGE(I,1,2)))
	      ELFM(I) = 0.5 * 
     &              (EL(INDEX_EDGE(I,1,1)) + EL(INDEX_EDGE(I,1,2)))
	   Endif
	Enddo 
C---------------------------------------------------------------------------------C
C              Define the wet and dry cells                                       c
c---------------------------------------------------------------------------------c
	Do I = 1, IJM
	   If(DC(I) .LE. 0.0) Then
	      CCM(I) = 0.0
	   Else
	      CCM(I) = 1.0
	   Endif
	Enddo
	Do I = 1, IJE
	   If(DS(I) .GT. 0.0 .AND. CFM(I) .EQ. 0.0) CFM(I) = 1.0
      Enddo
      
C	Do I = 1, IJE
C	   If(DS(I) .LE. 0.0) Then
C	      CFM(I) = 0.0
C	   Else
C	      CFM(I) = 1.0
C	   Endif
C	   If(INDEX_EDGE(I,1,1) .EQ. -999 .OR.
C     &      INDEX_EDGE(I,1,2) .EQ. -999) CFM(I) = -1.0
C	Enddo
      
	Do I = 1, IJM
	   If(CCM(I) .EQ. 0.0) Then
	      Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) 
     &			CFM(CELL_SIDE(I,J,1)) = 0.0
	      Enddo
	   Endif
	Enddo

      Do N = 1, EPTS
         ESAVE(N) = 0.0
	   Do K = 1, KBM
	      USAVE(N,K) = 0.0
	      VSAVE(N,K) = 0.0
	      WSAVE(N,K) = 0.0
	   Enddo
      Enddo

      Do N = 1, VPTS
         DZSAVE(N) = 0.0
	   UDSAVE(N) = 0.0
	   VDSAVE(N) = 0.0
         Do K = 1, KB
            UZSAVE(N,K) = 0.0
            VZSAVE(N,K) = 0.0
            SZSAVE(N,K) = 0.0
            TZSAVE(N,K) = 0.0
            VKM(N,K) = 0.0
            VKH(N,K) = 0.0
         Enddo
      Enddo
      Do K = 1, KB
         Do N = 1, FPTS
            CCFLUX(N,K) = 0.0
	   Enddo
      Enddo
      
      CONSPLT = .TRUE.
      CONSTSR = .TRUE.
c--------------------------------------------------------------------------C
      
      Return
      
      End
