C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (21/02/2023)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                                                                           #     
C               Ghost Cell Immersed Boundary Method moudle                  #                                           
C                                                                           #
C############################################################################
      Subroutine IBMALGC

      Include './Include/OCERM_INF'
      
      Real(kind=8) GCU(IJM,KBM),GCV(IJM,KBM),GCW(IJM,KBM)
      Real(kind=8) GCP(IJM,KBM),GCE(IJM)
      Real(kind=8) VEC_N(3),VEC_T(3),VEC_TEMP(3)
      Real(kind=8) UP_N(3),UP_T(3),UG_N(3),UG_T(3)
      Real(kind=8) UTEMP,VTEMP,WTEMP,PTEMP,ETEMP,YTEMP
      Real(kind=8) UTEMP2,VTEMP2,WTEMP2,PTEMP2,ETEMP2,YTEMP2
      Real(kind=8) ETAK,WM,DK,ZALPHA
      Real(kind=8) KAR,EWRLL,YIP,YGC,UN_MAG,UT_MAG,YPLUS,REST
      Integer NEIGHBORS,NUM
      
      
C===========================================================================C      
C         Initialize the arrays      
C===========================================================================C
      KAR = 0.4
      EWALL = 9.8
      NEIGHBORS = 28
      
      Do K = 1,KBM
      Do I = 1,IJM
          GCU(I,K) = 0
          GCV(I,K) = 0
          GCW(I,K) = 0
          GCP(I,K) = 0
          BFX(I,K) = 0
          BFY(I,K) = 0
          BFZ(I,K) = 0
          ZRC(I,K) = ZZ(K)*(HC(I)+ELF(I))+ELF(I)
      Enddo
      Enddo
      
      Do N = 1,IBNUM
      Do I = 1,5
          UGP(N,I) = 0
      Enddo
      Enddo

C===========================================================================C
C         Interpolation of image points
C===========================================================================C

      If (IBGCSCHEME.GT.0) Then 
          
      If (IBMDOMAIN.EQ.2) Then
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C           
      Do N = 1,NUMGH
      If (IDC(IGC(N,1),IGC(N,2)).EQ.2) Then  
          ETAK = 0
          Do J = 1,NEIGHBORS
          I = NBR(IGC(N,1),J)
          K = BGP(N,4)
          If (IDC(I,K).LT.0.9) Then
              DK = (BGP(N,1)-CXY(I,1))**2 + (BGP(N,2)-CXY(I,2))**2
              If (DK.LE.RITP**2) Then
                  If (DK.GT.1.0E-10) Then
                      ETAK = ETAK + 1/DK
                  Else
                      ETAK = ETAK + 1.0E10
                  Endif
              Endif
          Endif
          Enddo
          
          WM = 0
          UTEMP = 0
          VTEMP = 0
          WTEMP = 0
          PTEMP = 0
          ETEMP = 0
          
          Do J = 1,NEIGHBORS
          I = NBR(IGC(N,1),J)
          K = BGP(N,4)
          If (IDC(I,K).LT.0.9) Then
              DK = (BGP(N,1)-CXY(I,1))**2 + (BGP(N,2)-CXY(I,2))**2
              If (DK.LE.RITP**2) Then
                  If (DK.GT.1.0E-10) Then
                      WM = (1/ETAK) * (1/DK)
                  Else
                      WM = (1/ETAK) * 1.0E10
                  Endif
                  UTEMP = UTEMP + WM * U(I,K)
                  VTEMP = VTEMP + WM * V(I,K)
                  WTEMP = WTEMP + WM * QZ(I,K)
                  ETEMP = ETEMP + WM * ELF(I)
                  PTEMP = PTEMP + WM * PN(I,K)
              Endif
          Endif
          Enddo
          
          UGP(N,1) = UTEMP
          UGP(N,2) = VTEMP
          UGP(N,3) = WTEMP
          UGP(N,4) = ETEMP
          UGP(N,5) = PTEMP
      Else if (IDC(IGC(N,1),IGC(N,2)).GE.3) Then
          UGP(N,1) = U(BGP(N,1),BGP(N,4))
          UGP(N,2) = V(BGP(N,1),BGP(N,4))
          UGP(N,3) = QZ(BGP(N,1),BGP(N,4))
          UGP(N,4) = ELF(BGP(N,1))
          UGP(N,5) = PN(BGP(N,1),BGP(N,4))
      Endif
      Enddo
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      Endif
      
      If (IBMDOMAIN.EQ.3) Then
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Do N = 1,NUMGH
          K = BGP(N,4)
! Step1: K layer distance weight interpolation
          ETAK = 0
          Do J = 1,NEIGHBORS
          I = NBR(IGC(N,1),J)
          If (IDC(I,K).LT.0.9) Then
              DK = (BGP(N,1)-CXY(I,1))**2 + (BGP(N,2)-CXY(I,2))**2
              If (DK.LE.RITP**2) Then
                  If (DK.GT.1.0E-10) Then
                      ETAK = ETAK + 1/DK
                  Else
                      ETAK = ETAK + 1.0E10
                  Endif
              Endif
          Endif
          Enddo
          
          WM = 0
          UTEMP = 0
          VTEMP = 0
          WTEMP = 0
          PTEMP = 0
          ETEMP = 0
          
          Do J = 1,NEIGHBORS
          I = NBR(IGC(N,1),J)
          If (IDC(I,K).LT.0.9) Then
              DK = (BGP(N,1)-CXY(I,1))**2 + (BGP(N,2)-CXY(I,2))**2
              If (DK.LE.RITP**2) Then
                  If (DK.GT.1.0E-10) Then
                      WM = (1/ETAK) * (1/DK)
                  Else
                      WM = 0
                  Endif
                  UTEMP = UTEMP + WM * U(I,K)
                  VTEMP = VTEMP + WM * V(I,K)
                  WTEMP = WTEMP + WM * QZ(I,K)
                  ETEMP = ETEMP + WM * ELF(I)
                  PTEMP = PTEMP + WM * PN(I,K)
              Endif
          Endif
          Enddo
          
! Step2: K-1 layer distance weight interpolation          
          ETAK = 0
          Do J = 1,NEIGHBORS
          I = NBR(IGC(N,1),J)
          If (IDC(I,K-1).LT.0.9) Then
              DK = (BGP(N,1)-CXY(I,1))**2 + (BGP(N,2)-CXY(I,2))**2
              If (DK.LE.RITP**2) Then
                  If (DK.GT.1.0E-10) Then
                      ETAK = ETAK + 1/DK
                  Else
                      ETAK = ETAK + 1.0E10
                  Endif
              Endif
          Endif
          Enddo
          
          WM = 0
          UTEMP2 = 0
          VTEMP2 = 0
          WTEMP2 = 0
          PTEMP2 = 0
          ETEMP2 = 0
          
          Do J = 1,NEIGHBORS
          I = NBR(IGC(N,1),J)
          If (IDC(I,K-1).LT.0.9) Then
              DK = (BGP(N,1)-CXY(I,1))**2 + (BGP(N,2)-CXY(I,2))**2
              If (DK.LE.RITP**2) Then
                  If (DK.GT.1.0E-10) Then
                      WM = (1/ETAK) * (1/DK)
                  Else
                      WM = (1/ETAK) * 1.0E10
                  Endif
                  UTEMP2 = UTEMP2 + WM * U(I,K-1)
                  VTEMP2 = VTEMP2 + WM * V(I,K-1)
                  WTEMP2 = WTEMP2 + WM * QZ(I,K-1)
                  PTEMP2 = PTEMP2 + WM * PN(I,K-1)
              Endif
          Endif
          Enddo
          
! Step3: K & K-1 layer linear interpolation  
          I = IGC(N,1)
          ZALPHA = Abs(ZRC(I,K-1)-BGP(N,3)) / Abs(ZRC(I,K-1)-ZRC(I,K)) 
          UGP(N,1) = ZALPHA * UTEMP + (1 - ZALPHA) * UTEMP2
          UGP(N,2) = ZALPHA * VTEMP + (1 - ZALPHA) * VTEMP2
          UGP(N,3) = ZALPHA * WTEMP + (1 - ZALPHA) * WTEMP2
          UGP(N,4) = ETEMP
          UGP(N,5) = ZALPHA * PTEMP + (1 - ZALPHA) * PTEMP2
      Enddo
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      Endif
      
      Endif
      
C===========================================================================C      
C         Direct assignment method £¨Static boundary£©
C===========================================================================C
              
      If (IBGCSCHEME.EQ.0) Then 
      If (NSTEP.EQ.ISTART) Write(*,*) 'GCMODE is direct assignment.'
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      
      Do N = 1,NUMGH
          GCU(IGC(N,1),IGC(N,2)) = 0
          GCV(IGC(N,1),IGC(N,2)) = 0
          GCW(IGC(N,1),IGC(N,2)) = 0
      Enddo
      
      Do N = 1,NUMSL
          GCU(ISC(N,1),ISC(N,2)) = 0
          GCV(ISC(N,1),ISC(N,2)) = 0
          GCW(ISC(N,1),ISC(N,2)) = 0
      Enddo
      
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      Endif

C===========================================================================C
C         Distance Weight Linear interpolation
C===========================================================================C 
      
      If (IBGCSCHEME.EQ.1) Then
      If (NSTEP.EQ.ISTART) Then
          Write(*,*) 'GCMODE is linear weight interpolation.'
      Endif
      
      If (IBMDOMAIN.EQ.2) Then
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Do N = 1,NUMGH
      If (IDC(IGC(N,1),IGC(N,2)).EQ.2) Then  
          GCU(IGC(N,1),IGC(N,2)) = 0.5 * 
     &          (BXU(BGP(N,5),BGP(N,6),1) + UGP(N,1))
          GCV(IGC(N,1),IGC(N,2)) = 0.5 * 
     &          (BXU(BGP(N,5),BGP(N,6),2) + UGP(N,2))
          GCW(IGC(N,1),IGC(N,2)) = 0.5 * 
     &          (BXU(BGP(N,5),BGP(N,6),3) + UGP(N,3))
          GCE(IGC(N,1)) = UGP(N,4)
          GCP(IGC(N,1),IGC(N,2)) = UGP(N,5)
      Else if (IDC(IGC(N,1),IGC(N,2)).GE.3) Then
          If (BGP(N,10).LE.1.0E-10) Then
              ZALPHA = 1
          Else
              ZALPHA = 1 - Abs(BGP(N,3)-ZRC(IGC(N,1),IGC(N,2))) 
     &                 / BGP(N,10)
          Endif
          GCU(IGC(N,1),IGC(N,2)) = ZALPHA * UGP(N,1)
          GCV(IGC(N,1),IGC(N,2)) = ZALPHA * UGP(N,2)
          GCW(IGC(N,1),IGC(N,2)) = ZALPHA * UGP(N,3)
          GCE(IGC(N,1)) = UGP(N,4)
          GCP(IGC(N,1),IGC(N,2)) = UGP(N,5)
      Endif
      Enddo 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Endif
      
      If (IBMDOMAIN.EQ.3) Then
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Do N = 1,NUMGH
          GCU(IGC(N,1),IGC(N,2)) = 0.5 * 
     &          (BXU(BGP(N,5),BGP(N,6),1) + UGP(N,1))
          GCV(IGC(N,1),IGC(N,2)) = 0.5 * 
     &          (BXU(BGP(N,5),BGP(N,6),2) + UGP(N,2))
          GCW(IGC(N,1),IGC(N,2)) = 0.5 * 
     &          (BXU(BGP(N,5),BGP(N,6),3) + UGP(N,3))
          GCE(IGC(N,1)) = UGP(N,4)
          GCP(IGC(N,1),IGC(N,2)) = UGP(N,5)
      Enddo
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Endif

      Endif

C===========================================================================C
C         Wall function interpolation
C===========================================================================C 

      If (IBGCSCHEME.EQ.2) Then
      If (NSTEP.EQ.ISTART) Then
          Write(*,*) 'GCMODE is wall function interpolation.'
      Endif
      
      If (IBMDOMAIN.EQ.2) Then   
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Do N = 1,NUMGH  
          If (IDC(IGC(N,1),IGC(N,2)).EQ.2) Then      
              VEC_N(1) = BGP(N,7)
              VEC_N(2) = BGP(N,8)
              VEC_N(3) = 0
          Else if (IDC(IGC(N,1),IGC(N,2)).EQ.3) Then
              VEC_N(1) = 0
              VEC_N(2) = 0
              VEC_N(3) = -1
          Else if (IDC(IGC(N,1),IGC(N,2)).EQ.4) Then
              VEC_N(1) = 0
              VEC_N(2) = 0
              VEC_N(3) = 1
          Endif   
          
          YGC = BGP(N,10)
          YIP = BGP(N,10) * 2
          VEC_TEMP(1) = UGP(N,1)
          VEC_TEMP(2) = UGP(N,2)
          VEC_TEMP(3) = UGP(N,3)
          
          Call DOT_VEC(VEC_TEMP,VEC_N,UN_MAG)
          UP_N(1) = VEC_N(1) * UN_MAG
          UP_N(2) = VEC_N(2) * UN_MAG
          UP_N(3) = VEC_N(3) * UN_MAG
          UP_T(1) = VEC_TEMP(1) - UP_N(1)
          UP_T(2) = VEC_TEMP(2) - UP_N(2)
          UP_T(3) = 0
          UT_MAG = Sqrt(UP_T(1)**2 + UP_T(2)**2 + UP_T(3)**2)
              
          ! Calculate U normal of ghost cells
          UG_N(1) = UP_N(1) * YGC / (YIP + 1.0E-10)
          UG_N(2) = UP_N(2) * YGC / (YIP + 1.0E-10)
          UG_N(3) = UP_N(3) * YGC / (YIP + 1.0E-10)
      
          ! Calculate U tangent of ghost cells
          NUM = 0
          REST = 100
          YPLUS = 50
          Do while (REST.GT.1E-10) 
              YTEMP = YPLUS
              YPLUS = (YPLUS + KAR * YIP * UT_MAG / UMOL) / 
     &                (1 + log(EWALL * YPLUS))
              REST = Abs(YPLUS - YTEMP) 
              NUM = NUM + 1
              If (NUM.GT.1000) REST = 0
          End do
                    
          If (YPLUS.GE.15) Then
              UFR = YPLUS * UMOL / (YIP + 1.0E-10)
              
              If (UT_MAG.GT.0) Then
              UG_T(1) = (UT_MAG - (YIP - YGC) * UFR / (KAR * YPLUS)) * 
     &                (UP_T(1) / UT_MAG)
              UG_T(2) = (UT_MAG - (YIP - YGC) * UFR / (KAR * YPLUS)) * 
     &                (UP_T(2) / UT_MAG)
              UG_T(3) = (UT_MAG - (YIP - YGC) * UFR / (KAR * YPLUS)) * 
     &                (UP_T(3) / UT_MAG)
              Endif
          
              TKE(IGC(N,1),IGC(N,2)) = UFR**2/0.3
              TDISS(IGC(N,1),IGC(N,2)) = Sqrt((0.00008/(YIP**2))**2+
     &                (Sqrt(TKE(IGC(N,1),IGC(N,2)))/(0.219*YIP))**2)
      
              GCU(IGC(N,1),IGC(N,2)) = UG_N(1) + UG_T(1)
              GCV(IGC(N,1),IGC(N,2)) = UG_N(2) + UG_T(2)
              GCW(IGC(N,1),IGC(N,2)) = UG_N(3) + UG_T(3)
              GCE(IGC(N,1)) = UGP(N,4)
              GCP(IGC(N,1),IGC(N,2)) = UGP(N,5)
          Else
              If (IDC(IGC(N,1),IGC(N,2)).EQ.2) Then  
                  GCU(IGC(N,1),IGC(N,2)) = DWLI * UGP(N,1)
                  GCV(IGC(N,1),IGC(N,2)) = DWLI * UGP(N,2)
                  GCW(IGC(N,1),IGC(N,2)) = DWLI * UGP(N,3)
                  GCE(IGC(N,1)) = UGP(N,4)
                  GCP(IGC(N,1),IGC(N,2)) = UGP(N,5)
              Else if (IDC(IGC(N,1),IGC(N,2)).GE.3) Then
                  If (BGP(N,10).LE.1.0E-10) Then
                      ZALPHA = 1
                  Else
                      ZALPHA = 1 - Abs(BGP(N,3)-ZRC(IGC(N,1),IGC(N,2))) 
     &                        / BGP(N,10)
                  Endif
                  GCU(IGC(N,1),IGC(N,2)) = ZALPHA * UGP(N,1)
                  GCV(IGC(N,1),IGC(N,2)) = ZALPHA * UGP(N,2)
                  GCW(IGC(N,1),IGC(N,2)) = ZALPHA * UGP(N,3)
                  GCE(IGC(N,1)) = UGP(N,4)
                  GCP(IGC(N,1),IGC(N,2)) = UGP(N,5)
              Endif
          Endif
      Enddo
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Endif
      
      If (IBMDOMAIN.EQ.3) Then  
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      Do N = 1,NUMGH
          VEC_N(1) = BGP(N,7)
          VEC_N(2) = BGP(N,8)
          VEC_N(3) = BGP(N,9)
          YGC = BGP(N,10)
          YIP = BGP(N,10) * 2
               
          VEC_TEMP(1) = UGP(N,1)
          VEC_TEMP(2) = UGP(N,2)
          VEC_TEMP(3) = UGP(N,3)
          Call DOT_VEC(VEC_TEMP,VEC_N,UN_MAG)
          UP_N(1) = VEC_N(1) * UN_MAG
          UP_N(2) = VEC_N(2) * UN_MAG
          UP_N(3) = VEC_N(3) * UN_MAG
          UP_T(1) = UTEMP - UP_N(1)
          UP_T(2) = VTEMP - UP_N(2)
          UP_T(3) = WTEMP - UP_N(3)
          UT_MAG = Sqrt(UP_T(1)**2 + UP_T(2)**2 + UP_T(3)**2)
      
          ! Calculate U normal of ghost cells
          UG_N(1) = UP_N(1) * YGC / YIP
          UG_N(2) = UP_N(2) * YGC / YIP
          UG_N(3) = UP_N(3) * YGC / YIP
      
          ! Calculate U tangent of ghost cells
          NUM = 0
          REST = 100
          YPLUS = 20
          Do while (REST.GT.1E-10) 
              YTEMP = YPLUS
              YPLUS = (YPLUS + KAR * YIP * UT_MAG / UMOL) / 
     &                (1 + log(EWALL * YPLUS))
              REST = Abs(YPLUS - YTEMP) 
              NUM = NUM + 1
              If (NUM.GT.1000) REST = 0
          End do
          UFR = YPLUS * UMOL / YIP
      
          UG_T(1) = (UT_MAG - (YIP - YGC) * UFR / (KAR * YPLUS)) * 
     &              (UP_T(1) / UT_MAG)
          UG_T(2) = (UT_MAG - (YIP - YGC) * UFR / (KAR * YPLUS)) * 
     &              (UP_T(2) / UT_MAG)
          UG_T(3) = (UT_MAG - (YIP - YGC) * UFR / (KAR * YPLUS)) * 
     &              (UP_T(3) / UT_MAG)
      
          TKE(IGC(N,1),IGC(N,2)) = UFR**2/0.3
          TDISS(IGC(N,1),IGC(N,2)) = Sqrt((0.00008/(YIP**2))**2+
     &             (Sqrt(TKE(IGC(N,1),IGC(N,2)))/(0.219*YIP))**2)

          GCU(IGC(N,1),IGC(N,2)) = UG_N(1) + UG_T(1)
          GCV(IGC(N,1),IGC(N,2)) = UG_N(2) + UG_T(2)
          GCW(IGC(N,1),IGC(N,2)) = UG_N(3) + UG_T(3)
          GCE(IGC(N,1)) = UGP(N,4)
          GCP(IGC(N,1),IGC(N,2)) = UGP(N,5)
      Enddo
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Endif
      
      Endif

C===========================================================================C      
C         Update value to fresh cells
C===========================================================================C

      If (IBMSCHEME.GE.1) Then
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Do K = 1,KBM
      Do I = 1,IJM
      If (IFC(I,K).GT.0.9) Then
          ETAK = 1 / BGP(N,10)
          Do J = 1,CELL_POLYGEN(I)
          If (IDC(CELL_SIDE(I,J,2),K).LT.0.9) Then
              DK = (CXY(CELL_SIDE(I,J,2),1)-CXY(I,1))**2 +
     &             (CXY(CELL_SIDE(I,J,2),2)-CXY(I,2))**2
              ETAK = ETAK + 1/DK
          Endif
          Enddo
          
          WM = 0
          UTEMP = BXU(BGP(N,5),BGP(N,6),1) * (1/ETAK) * (1/BGP(N,10))
          VTEMP = BXU(BGP(N,5),BGP(N,6),2) * (1/ETAK) * (1/BGP(N,10))
          WTEMP = BXU(BGP(N,5),BGP(N,6),3) * (1/ETAK) * (1/BGP(N,10))
          PTEMP = 0
          Do J = 1,CELL_POLYGEN(I)
          If (IDC(CELL_SIDE(I,J,2),K).LT.0.9) Then
              DK = (CXY(CELL_SIDE(I,J,2),1)-CXY(I,1))**2 +
     &             (CXY(CELL_SIDE(I,J,2),2)-CXY(I,2))**2
              WM = (1/ETAK) * (1/DK)
              UTEMP = UTEMP + WM * U(CELL_SIDE(I,J,2),K)
              VTEMP = VTEMP + WM * V(CELL_SIDE(I,J,2),K)
              WTEMP = WTEMP + WM * QZ(CELL_SIDE(I,J,2),K)
              PTEMP = PTEMP + WM * PN(CELL_SIDE(I,J,2),K)
          Endif
          Enddo
          
          U(I,K) = UTEMP
          V(I,K) = VTEMP
          QZ(I,K) = WTEMP
          PN(I,K) = PTEMP
      Endif
      Enddo
      Enddo
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Endif
C===========================================================================C      
C         Assign value to solid cells and ghost cells
C===========================================================================C
            
      Do K = 1, KBM
      Do I = 1, IJM
          If (IDC(I,K).GT.0.9) Then
              BFX(I,K) = (GCU(I,K) - U(I,K)) / DTI
              BFY(I,K) = (GCV(I,K) - V(I,K)) / DTI
              BFZ(I,K) = (GCW(I,K) - QZ(I,K)) / DTI           
          Endif
      Enddo
      Enddo
      
      Do K = 1, KBM
      Do I = 1, IJM
          If (IDC(I,K).GT.0.9) Then
              U(I,K) = GCU(I,K)
              V(I,K) = GCV(I,K)
              QZ(I,K) = GCW(I,K)
              If(IDC(I,K).EQ.1) Then
                  TKE(I,K) = 0
                  TDISS(I,K) = 100000
              Endif  
          Endif
      Enddo
      Enddo
      
      If (IBGCMODE.EQ.2) Then
          Do K = 1, KBM
          Do I = 1, IJM
              If (IDC(I,K).GT.0.9) Then
                  PN(I,K) = GCP(I,K)  
              Endif
          Enddo
          Enddo
      Endif
      
1000  Format (3F8.4)
1100  Format (4F10.6)
      
      
      End Subroutine