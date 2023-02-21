C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (21/02/2023)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c                                                                           #
c                 Collision Detection of Rigid Bodies                       #                                           
c                                                                           #
c############################################################################

	Subroutine IBMCLBD(BODY,RGROT,RGI1,BODY2)

      Include './Include/OCERM_INF'
	  
      Integer BODY,BODY2
      Real(Kind=8) RGROT(3,3),RGI1(3,3)
      Real(Kind=8) MAT_RRI(3,3),MAT_K(3,3),MAT_K1(3,3)
      Real(kind=8) CRP(3),BDT(3),PLANE(3,3),BXY3(10000,3,3)
      Real(Kind=8) VEC_N(3),VEC_J(3),RI(3),RRI(3),WRI(3)
      Real(Kind=8) UI(3),UI1(3),UI_N(3),UI_T(3),UI_SUM(3),RI_SUM(3)
      Real(Kind=8) MAT_TEMP(3,3),VEC_TEMP(3),VEC_TEMP2(3)
      Real(Kind=8) XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN,FL,FR,GK,ZIB
      Real(Kind=8) PA,PB,PC,PD,PDIS,PDMIN,TLINE,FRT,N_MAG,N_SCA
      Integer CLP(10000),INT,NCL,NTC,NUM_CL
      Logical is_coincide
      
      
      Do I = 1,10000
      Do K = 1,3
          BXY3(I,K,1) = 0
          BXY3(I,K,2) = 0
          BXY3(I,K,3) = 0
      Enddo
      CLP(I) = 0
      Enddo
      
      Do I = 1,3
          RI_SUM(I) = 0
          UI_SUM(I) = 0
      Enddo
      NUM_CL = 0
      NCL = 0
      
      NTC = 0     ! NTC is used to count triangle element
      
      Do I = 1,IBN(BODY2)/3
          NTC = NTC+1
              
          BXY3(NTC,1,1) = BXY(BODY2,3*I-2,1)
          BXY3(NTC,1,2) = BXY(BODY2,3*I-2,2)
          BXY3(NTC,1,3) = BXY(BODY2,3*I-2,3)

          BXY3(NTC,2,1) = BXY(BODY2,3*I-1,1)
          BXY3(NTC,2,2) = BXY(BODY2,3*I-1,2)
          BXY3(NTC,2,3) = BXY(BODY2,3*I-1,3)

          BXY3(NTC,3,1) = BXY(BODY2,3*I,1)
          BXY3(NTC,3,2) = BXY(BODY2,3*I,2)
          BXY3(NTC,3,3) = BXY(BODY2,3*I,3)
      Enddo

      RGJ(BODY2,1) = 0
          RGJ(BODY2,2) = 0
          RGJ(BODY2,3) = 0
          RGJ(BODY2+2,1) = 0
          RGJ(BODY2+2,2) = 0
          RGJ(BODY2+2,3) = 0
          
C===========================================================================c     
C         Set Bounding Box     
C===========================================================================c

      XMAX = BXY(BODY2,1,1)
      XMIN = BXY(BODY2,1,1)
      YMAX = BXY(BODY2,1,2)
      YMIN = BXY(BODY2,1,2)
      ZMAX = BXY(BODY2,1,3)
      ZMIN = BXY(BODY2,1,3)

      Do I = 1,IBN(BODY2)
          If (BXY(BODY2,I,1).GE.XMAX) Then
              XMAX = BXY(BODY2,I,1)
          Endif
          If (BXY(BODY2,I,1).LE.XMIN) Then
              XMIN = BXY(BODY2,I,1)
          Endif
          If (BXY(BODY2,I,2).GE.YMAX) Then
              YMAX = BXY(BODY2,I,2)
          Endif
          If (BXY(BODY2,I,2).LE.YMIN) Then
              YMIN = BXY(BODY2,I,2)
          Endif
          If (BXY(BODY2,I,3).GE.ZMAX) Then
              ZMAX = BXY(BODY2,I,3)
          Endif
          If (BXY(BODY2,I,3).LE.ZMIN) Then
              ZMIN = BXY(BODY2,I,3)
          Endif
      Enddo 
      
C===========================================================================c     
C         Ray-casting Method       
C===========================================================================c
      
      Do 300 N = 1,IBN(BODY)
          
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      BDT(1) = BXY(BODY,N,1)
      BDT(2) = BXY(BODY,N,2)
      BDT(3) = BXY(BODY,N,3) 
               
      If ((BDT(1).LE.XMAX).AND.(BDT(1).GE.XMIN).AND. 
     &    (BDT(2).LE.YMAX).AND.(BDT(2).GE.YMIN).AND.
     &    (BDT(3).LE.ZMAX).AND.(BDT(3).GE.ZMIN)) Then 
      INT = 0
      Do 301 I = 1,NTC
      GK = (BXY3(I,2,1)-BXY3(I,1,1))*(BXY3(I,3,2)-BXY3(I,1,2))
     &    -(BXY3(I,2,2)-BXY3(I,1,2))*(BXY3(I,3,1)-BXY3(I,1,1))
      If (GK.NE.0) Then
          is_coincide = (((BDT(1).EQ.BXY3(I,1,1)).AND.
     &                (BDT(2).EQ.BXY3(I,1,2))).OR.
     &                ((BDT(1).EQ.BXY3(I,2,1)).AND.
     &                (BDT(2).EQ.BXY3(I,2,2))).OR.
     &                ((BDT(1).EQ.BXY3(I,3,1)).AND.
     &                (BDT(2).EQ.BXY3(I,3,2))))
          If (is_coincide) Then
              INT = 999
              Exit
          Endif

          CRP(1) = (BXY3(I,1,1)-BDT(1))*(BXY3(I,2,2)-BDT(2))
     &            -(BXY3(I,1,2)-BDT(2))*(BXY3(I,2,1)-BDT(1))
          CRP(2) = (BXY3(I,2,1)-BDT(1))*(BXY3(I,3,2)-BDT(2))
     &            -(BXY3(I,2,2)-BDT(2))*(BXY3(I,3,1)-BDT(1))
          CRP(3) = (BXY3(I,3,1)-BDT(1))*(BXY3(I,1,2)-BDT(2))
     &            -(BXY3(I,3,2)-BDT(2))*(BXY3(I,1,1)-BDT(1))

          If (((CRP(1).GE.0).AND.(CRP(2).GE.0).AND.(CRP(3).GE.0)).OR.
     &        ((CRP(1).LT.0).AND.(CRP(2).LT.0).AND.(CRP(3).LT.0))) Then

          FL=((BXY3(I,2,1)-BXY3(I,1,1))*(BXY3(I,3,2)-BXY3(I,1,2))-
     &        (BXY3(I,2,2)-BXY3(I,1,2))*(BXY3(I,3,1)-BXY3(I,1,1))) 
          FR =(((BDT(2)-BXY3(I,1,2))*(BXY3(I,2,1)-BXY3(I,1,1))*
     &        (BXY3(I,3,3)-BXY3(I,1,3)))+
     &        ((BDT(1)-BXY3(I,1,1))*(BXY3(I,2,3)-BXY3(I,1,3))*
     &        (BXY3(I,3,2)-BXY3(I,1,2)))-
     &        ((BDT(1)-BXY3(I,1,1))*(BXY3(I,2,2)-BXY3(I,1,2))*
     &        (BXY3(I,3,3)-BXY3(I,1,3)))-
     &        ((BDT(2)-BXY3(I,1,2))*(BXY3(I,2,3)-BXY3(I,1,3))*
     &        (BXY3(I,3,1)-BXY3(I,1,1))))
          ZIB = FR/FL+BXY3(I,1,3)
          
          If (BDT(3).LT.ZIB) Then
              INT = INT + 1
          Else If (BDT(3).EQ.ZIB) Then
              INT = 1
              Exit
          Endif
          Endif
      Endif 
                   
301   Continue
      
      If (mod(INT,2) .GT. 0.0) Then
          NCL = NCL + 1
          CLP(NCL) = N
      Endif
      
      Endif 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
       
300   Continue
      
C===========================================================================c     
C         Impuslse Method  
C===========================================================================c 
      
      If (NCL.EQ.0) Return
      
      Do M = 1,NCL
          N = CLP(M)
          PDMIN = 10000
          Do I = 1,NTC
		    PDIS = (BXY(BODY,N,1) - (BXY3(I,1,1) +
     &               BXY3(I,2,1) + BXY3(I,3,1))/3)**2 +
     &               (BXY(BODY,N,2) - (BXY3(I,1,2) +
     &               BXY3(I,2,2) + BXY3(I,3,2))/3)**2 +
     &               (BXY(BODY,N,3) - (BXY3(I,1,3) +
     &               BXY3(I,2,3) + BXY3(I,3,3))/3)**2
		    If (PDIS.LE.PDMIN) Then
			    PDMIN = PDIS
                  PLANE(1,1) = BXY3(I,1,1)
                  PLANE(1,2) = BXY3(I,1,2)
                  PLANE(1,3) = BXY3(I,1,3)
                  PLANE(2,1) = BXY3(I,2,1)
                  PLANE(2,2) = BXY3(I,2,2)
                  PLANE(2,3) = BXY3(I,2,3)
                  PLANE(3,1) = BXY3(I,3,1)
                  PLANE(3,2) = BXY3(I,3,2)
                  PLANE(3,3) = BXY3(I,3,3)
              Endif
          Enddo 
          PA = (PLANE(2,2)-PLANE(1,2)) * (PLANE(3,3)-PLANE(1,3))
     &        -(PLANE(2,3)-PLANE(1,3)) * (PLANE(3,2)-PLANE(1,2))
          PB = (PLANE(2,3)-PLANE(1,3)) * (PLANE(3,1)-PLANE(1,1))
     &        -(PLANE(2,1)-PLANE(1,1)) * (PLANE(3,3)-PLANE(1,3))
          PC = (PLANE(2,1)-PLANE(1,1)) * (PLANE(3,2)-PLANE(1,2))
     &        -(PLANE(2,2)-PLANE(1,2)) * (PLANE(3,1)-PLANE(1,1))
          PD = -1 * (PA*PLANE(1,1)+PB*PLANE(1,2)+PC*PLANE(1,3))
          TLINE = -1 * (PA*BXY(BODY,N,1) + PB*BXY(BODY,N,2) +
     &            PC*BXY(BODY,N,3) + PD) / (PA**2 + PB**2 + PC**2)  
          VEC_N(1) = PA * TLINE
          VEC_N(2) = PB * TLINE
          VEC_N(3) = PC * TLINE
          N_MAG = Sqrt(VEC_N(1)**2 + VEC_N(2)**2 + VEC_N(3)**2) 
     &            + 0.0000001
          VEC_N(1) = VEC_N(1) / N_MAG
          VEC_N(2) = VEC_N(2) / N_MAG
          VEC_N(3) = VEC_N(3) / N_MAG
          
          RI(1) = RI0(BODY,N,1)
          RI(2) = RI0(BODY,N,2)
          RI(3) = RI0(BODY,N,3)
          Call MUL_VEC(RGROT,RI,RRI)
          VEC_TEMP(1) = RGW(BODY,1)
          VEC_TEMP(2) = RGW(BODY,2)
          VEC_TEMP(3) = RGW(BODY,3)
          Call CRP_VEC(VEC_TEMP,RRI,WRI)
          UI(1) = RGU(BODY,1) + WRI(1)
          UI(2) = RGU(BODY,2) + WRI(2)
          UI(3) = RGU(BODY,3) + WRI(3)
          Call DOT_VEC(UI,VEC_N,N_SCA)
          If (N_SCA.LE.0) Then
              RI_SUM(1) = RI_SUM(1) + RRI(1)
              RI_SUM(2) = RI_SUM(2) + RRI(2)
              RI_SUM(3) = RI_SUM(3) + RRI(3)
              UI_SUM(1) = UI_SUM(1) + UI(1)
              UI_SUM(2) = UI_SUM(2) + UI(2)
              UI_SUM(3) = UI_SUM(3) + UI(3)
              NUM_CL = NUM_CL + 1
          Endif
      Enddo
      
      If (NUM_CL.GT.0) Then
          ! Compute relative position and velocity of collision point
          RRI(1) = RI_SUM(1) / NUM_CL
          RRI(2) = RI_SUM(2) / NUM_CL
          RRI(3) = RI_SUM(3) / NUM_CL
          UI(1) = UI_SUM(1) / NUM_CL
          UI(2) = UI_SUM(2) / NUM_CL
          UI(3) = UI_SUM(3) / NUM_CL
          
          ! Compute the normal vector of collision plane
          RI(1) = RGX(BODY,1)
          RI(2) = RGX(BODY,2)
          RI(3) = RGX(BODY,3)
          PDMIN = 10000
          Do I = 1,NTC
		    PDIS = (RI(1) - (BXY3(I,1,1) + 
     &               BXY3(I,2,1) + BXY3(I,3,1))/3)**2 +
     &               (RI(2) - (BXY3(I,1,2) +
     &               BXY3(I,2,2) + BXY3(I,3,2))/3)**2 +
     &               (RI(3) - (BXY3(I,1,3) +
     &               BXY3(I,2,3) + BXY3(I,3,3))/3)**2
		    If (PDIS.LE.PDMIN) Then
			    PDMIN = PDIS
                  PLANE(1,1) = BXY3(I,1,1)
                  PLANE(1,2) = BXY3(I,1,2)
                  PLANE(1,3) = BXY3(I,1,3)
                  PLANE(2,1) = BXY3(I,2,1)
                  PLANE(2,2) = BXY3(I,2,2)
                  PLANE(2,3) = BXY3(I,2,3)
                  PLANE(3,1) = BXY3(I,3,1)
                  PLANE(3,2) = BXY3(I,3,2)
                  PLANE(3,3) = BXY3(I,3,3)
              Endif
          Enddo 
          PA = (PLANE(2,2)-PLANE(1,2)) * (PLANE(3,3)-PLANE(1,3))
     &        -(PLANE(2,3)-PLANE(1,3)) * (PLANE(3,2)-PLANE(1,2))
          PB = (PLANE(2,3)-PLANE(1,3)) * (PLANE(3,1)-PLANE(1,1))
     &        -(PLANE(2,1)-PLANE(1,1)) * (PLANE(3,3)-PLANE(1,3))
          PC = (PLANE(2,1)-PLANE(1,1)) * (PLANE(3,2)-PLANE(1,2))
     &        -(PLANE(2,2)-PLANE(1,2)) * (PLANE(3,1)-PLANE(1,1))
          PD = -1 * (PA*PLANE(1,1)+PB*PLANE(1,2)+PC*PLANE(1,3))
          TLINE = -1 * (PA*RI(1) + PB*RI(2) + PC*RI(3) + PD) /
     &            (PA**2 + PB**2 + PC**2) 
          VEC_N(1) = PA * TLINE
          VEC_N(2) = PB * TLINE
          VEC_N(3) = PC * TLINE
          N_MAG = Sqrt(VEC_N(1)**2 + VEC_N(2)**2 + VEC_N(3)**2) 
     &            + 0.0000001
          VEC_N(1) = VEC_N(1) / N_MAG
          VEC_N(2) = VEC_N(2) / N_MAG
          VEC_N(3) = VEC_N(3) / N_MAG
          
          ! Compute nomral and tangent velocity of collision point
          Call DOT_VEC(UI,VEC_N,N_SCA)
          UI_N(1) = VEC_N(1) * N_SCA
          UI_N(2) = VEC_N(2) * N_SCA
          UI_N(3) = VEC_N(3) * N_SCA
          UI_T(1) = UI(1) - UI_N(1)
          UI_T(2) = UI(2) - UI_N(2)
          UI_T(3) = UI(3) - UI_N(3)

          ! Set new relative velocity of collision point
          FRT = 1 - 0.8 * (1 + FRN(BODY)) * 
     &          (Sqrt(UI_N(1)**2 + UI_N(2)**2 + UI_N(3)**2)) / 
     &          (Sqrt(UI_T(1)**2 + UI_T(2)**2 + UI_T(3)**2)
     &           + 0.00000001)
          If (FRT.LE.0) FRT = 0
          UI1(1) = FRT * UI_T(1) - FRN(BODY) * UI_N(1)
          UI1(2) = FRT * UI_T(2) - FRN(BODY) * UI_N(2)
          UI1(3) = FRT * UI_T(3) - FRN(BODY) * UI_N(3)
          
          ! Set star matrix of collision point
          MAT_RRI(1,1) = 0
          MAT_RRI(1,2) = -1 * RRI(3)
          MAT_RRI(1,3) = RRI(2)
          MAT_RRI(2,1) = RRI(3)
          MAT_RRI(2,2) = 0
          MAT_RRI(2,3) = -1 * RRI(1)
          MAT_RRI(3,1) = -1 * RRI(2)
          MAT_RRI(3,2) = RRI(1)
          MAT_RRI(3,3) = 0
          
          ! Compute matrix K
          MAT_TEMP = MATMUL(MAT_RRI,RGI1)
          MAT_K = MATMUL(MAT_TEMP,MAT_RRI)
          
          Do I = 1,3
          Do J = 1,3
              If (I.EQ.J) Then
                  MAT_K(I,J) = 1 / RGM(BODY) - MAT_K(I,J)
              Else
                  MAT_K(I,J) = 0 - MAT_K(I,J)
              Endif
          Enddo
          Enddo
          
          ! Compute impulse J
          Call MAT_INV(MAT_K,MAT_K1)
          VEC_TEMP(1) = UI1(1) - UI(1) 
          VEC_TEMP(2) = UI1(2) - UI(2)
          VEC_TEMP(3) = UI1(3) - UI(3)
          Call MUL_VEC(MAT_K1,VEC_TEMP,VEC_J) 
          RGJ(BODY2,1) = -1 * VEC_J(1)
          RGJ(BODY2,2) = -1 * VEC_J(2)
          RGJ(BODY2,3) = -1 * VEC_J(3)
          RGJ(BODY2+2,1) = RGX(BODY,1) + RRI(1)
          RGJ(BODY2+2,2) = RGX(BODY,2) + RRI(2)
          RGJ(BODY2+2,3) = RGX(BODY,3) + RRI(3)
          ! Update u and w of BODY
          RGU(BODY,1) = RGU(BODY,1) + 1 / RGM(BODY) * VEC_J(1)
          RGU(BODY,2) = RGU(BODY,2) + 1 / RGM(BODY) * VEC_J(2)
          RGU(BODY,3) = RGU(BODY,3) + 1 / RGM(BODY) * VEC_J(3)
          Call CRP_VEC(RRI,VEC_J,VEC_TEMP)
          Call MUL_VEC(RGI1,VEC_TEMP,VEC_TEMP2) 
          RGW(BODY,1) = RGW(BODY,1) + VEC_TEMP2(1)
          RGW(BODY,2) = RGW(BODY,2) + VEC_TEMP2(2)
          RGW(BODY,3) = RGW(BODY,3) + VEC_TEMP2(3)
      Endif

      End Subroutine