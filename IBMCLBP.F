C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (21/02/2023)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c                                                                           #
c                 Collision Detection of PLANE Boundary                     #                                           
c                                                                           #
c############################################################################
      Subroutine IBMCLBP(BODY,RGROT,RGI1,PLANE)

      Include './Include/OCERM_INF'
      
      Integer BODY
      Real(Kind=8) RGROT(3,3),RGI1(3,3),PLANE(4)
      Real(Kind=8) MAT_RRI(3,3),MAT_K(3,3),MAT_K1(3,3)
      Real(Kind=8) VEC_N(3),VEC_J(3),RI_SUM(3),UI_SUM(3)
      Real(Kind=8) RI(3),RRI(3),WRI(3),UI(3),UI1(3),UI_N(3),UI_T(3)
      Real(Kind=8) FRT,PL_MAG,N_SCA,BIAS
      Real(Kind=8) MAT_TEMP(3,3),VEC_TEMP(3),VEC_TEMP2(3)
      Integer CLP(10000),IND_CL,INT,NCL,NTC,NUM_CL
      
      
      Do I = 1,3
      Do J = 1,3
          MAT_K(I,J) = 0
          MAT_K1(I,J) = 0
      Enddo
      Enddo
      
      Do I = 1,10000
          CLP(I) = 0
      Enddo
      
      Do I = 1,3
          VEC_J(I) = 0
          RI_SUM(I) = 0
          UI_SUM(I) = 0
      Enddo
      
      VEC_N(1) = PLANE(1)
      VEC_N(2) = PLANE(2)
      VEC_N(3) = PLANE(3)
      
      BIAS = 0.000001
C===========================================================================c     
C         Impuslse Method  
C===========================================================================c 
      
      PL_MAG = PLANE(1) + PLANE(2) + PLANE(3)
      If (PL_MAG.GT.0) Then
          If (PLANE(1).NE.0) IND_CL = 1 
          If (PLANE(2).NE.0) IND_CL = 2 
          If (PLANE(3).NE.0) IND_CL = 3 
          NCL = 0
          Do N = 1,IBN(BODY)  
          ! Collision points detection
              If (BXY(BODY,N,IND_CL).LT.PLANE(4)) Then 
                  NCL = NCL + 1
                  CLP(NCL) = N
              Endif
          Enddo
      Else
          If (PLANE(1).NE.0) IND_CL = 1 
          If (PLANE(2).NE.0) IND_CL = 2 
          If (PLANE(3).NE.0) IND_CL = 3 
          NCL = 0
          Do N = 1,IBN(BODY)  
          ! Collision points detection
              If (BXY(BODY,N,IND_CL).GT.PLANE(4)) Then
                  NCL = NCL + 1
                  CLP(NCL) = N
              Endif
          Enddo
      Endif
          
c      If (NCL.EQ.0) Return
              
      Do I = 1,3
          RI_SUM(I) = 0
          UI_SUM(I) = 0
      Enddo
      NUM_CL = 0
      
      Do M = 1,NCL    
          RI(1) = RI0(BODY,CLP(M),1)
          RI(2) = RI0(BODY,CLP(M),2)
          RI(3) = RI0(BODY,CLP(M),3)
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
      
      Do I = 1,3
          RRI(I) = 0
          WRI(I) = 0
      Enddo
      
      If (NUM_CL.GT.0) Then
          ! Compute relative position and velocity of collision point
          RRI(1) = RI_SUM(1) / NUM_CL
          RRI(2) = RI_SUM(2) / NUM_CL
          RRI(3) = RI_SUM(3) / NUM_CL 
          UI(1) = UI_SUM(1) / NUM_CL
          UI(2) = UI_SUM(2) / NUM_CL
          UI(3) = UI_SUM(3) / NUM_CL
          
          ! Compute nomral and tangent velocity of collision point
          Call DOT_VEC(UI,VEC_N,N_SCA)
          UI_N(1) = VEC_N(1) * N_SCA
          UI_N(2) = VEC_N(2) * N_SCA
          UI_N(3) = VEC_N(3) * N_SCA
          UI_T(1) = UI(1) - UI_N(1)
          UI_T(2) = UI(2) - UI_N(2)
          UI_T(3) = UI(3) - UI_N(3)
          If (Abs(UI_N(IND_CL)).LE.BIAS) FRN(BODY) = 0
          
          ! Set new relative velocity of collision point
          FRT = 1 - 0.8 * (1 + FRN(BODY)) * 
     &          (Sqrt(UI_N(1)**2 + UI_N(2)**2 + UI_N(3)**2)) / 
     &          (Sqrt(UI_T(1)**2 + UI_T(2)**2 + UI_T(3)**2)
     &           + 0.00000001)
          If ((FRT.LE.0).OR.(FRN(BODY).EQ.0)) FRT = 0
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
       
          ! Update u and w
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