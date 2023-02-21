C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (21/02/2023)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                                                                           #
C        Detecte differnet kind of cells and interpolation points.          #                                           
C                                                                           #
C############################################################################
      Subroutine IBMALIDC
      
      Include './Include/OCERM_INF'
      
      Real(kind=8) BXY2(10000,2,2),BXY3(10000,3,3)
      Real(kind=8) PLANE(3,3),LINE(2,2),CRP(3),Z2D(2)
      Double precision XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN
      Double precision FL,FR,GK,ZIB,XIM,YIM,ZIM,ZNE,RRR
      Double precision PA,PB,PC,PD,PDIS,PDMIN,TLINE,N_MAG,N_SCA
      Double precision IPX,IPY,IPZ,IPK,IPDIS,IPNV(3)
      Integer S1,N_SL,N_GH,NBMIN,NPMIN,P1,P2,P3
      Logical is_coincide
      

C===========================================================================C      
C         Initialize the arrays      
C===========================================================================C

      NUMSL = 0
      NUMGH = 0
      Z2D(1) = -0.13
      Z2D(2) = -0.06
      
      Do K = 1,KBM
      Do I = 1,IJM
          IFC(I,K) = IDC(I,K)
          IDC(I,K) = 0
          IDE(I,K) = 0
          ZRC(I,K) = ZZ(K)*(HC(I)+ELF(I))+ELF(I)
      Enddo
      Enddo 
      
      Do N = 1,IBNUM
          IGC(N,1) = 0
          IGC(N,2) = 0
          Do J = 1,9
              BGP(N,J) = 0
          Enddo
      Enddo  
      
      Do N = 1,IBNUM*2
          ISC(N,1) = 0
          ISC(N,2) = 0
      Enddo
      
      Do I = 1,3
          CRP(I) = 0
          IPNV(I) = 0
      Enddo
      
C===========================================================================C      
C         2D Immersed boundary detection      
C===========================================================================C   
      
      If (IBMDOMAIN.EQ.2) Then 
          
      Do 200 NB = 1,NUMIB
        
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      Do I = 1,10000
      Do K = 1,2
          BXY2(I,K,1) = 0
          BXY2(I,K,2) = 0
      Enddo
      Enddo
          
      Do I = 1, NUME2(NB)
          BXY2(I,1,1) = BXY(NB,IBE2(NB,I,1),1)
          BXY2(I,1,2) = BXY(NB,IBE2(NB,I,1),2)
          BXY2(I,2,1) = BXY(NB,IBE2(NB,I,2),1)
          BXY2(I,2,2) = BXY(NB,IBE2(NB,I,2),2)
      Enddo
      
! Get XMAX,XMIN,YMAX,YMIN 
      XMAX = BXY(NB,1,1)
      XMIN = BXY(NB,1,1)
      YMAX = BXY(NB,1,2)
      YMIN = BXY(NB,1,2)
     
      Do I = 1,IBN(NB)
          If (BXY(NB,I,1).GE.XMAX) Then
              XMAX = BXY(NB,I,1)
          Endif
          If (BXY(NB,I,1).LE.XMIN) Then
              XMIN = BXY(NB,I,1)
          Endif
          If (BXY(NB,I,2).GE.YMAX) Then
              YMAX = BXY(NB,I,2)
          Endif
          If (BXY(NB,I,2).LE.YMIN) Then
              YMIN = BXY(NB,I,2)
          Endif
      Enddo 

! Solid cells detection   
      Do 201 I = 1,IJM 
C---------------------------------------------------------------------------C
C             Ray-tracing Method                                            C
C---------------------------------------------------------------------------C            
      If ((CXY(I,1).LT.XMAX).AND.(CXY(I,1).GT.XMIN).AND. 
     &    (CXY(I,2).LT.YMAX).AND.(CXY(I,2).GT.YMIN)) Then
      S1 = 0   
      Do 202 N = 1,NUME2(NB)
          CRP(1) = (CXY(I,1)-BXY2(N,1,1)) * (BXY2(N,2,2)-CXY(I,2))
     &            -(BXY2(N,2,1)-CXY(I,1)) * (CXY(I,2)-BXY2(N,1,2))
                           
          If (((BXY2(N,1,1).LE.CXY(I,1)).AND.
     &        (BXY2(N,2,1).GT.CXY(I,1)).AND.(CRP(1).GE.0))
     &        .OR.
     &        ((BXY2(N,1,1).GE.CXY(I,1)).AND.
     &        (BXY2(N,2,1).LT.CXY(I,1)).AND.(CRP(1).LE.0))) 
     &    Then               
              S1 = S1 + 1
              If (CRP(1) .EQ. 0) Then    
                  S1 = 1
                  Exit
              Endif
          Endif  
202   Continue 
      
      Do K = 1,KBM
      If ((ZRC(I,K).GE.Z2D(1)).AND.(ZRC(I,K).LE.Z2D(2))) Then
          If (mod(S1,2) .GT. 0) Then
              IDC(I,K) = 1     
          Endif
      Endif
      Enddo
      
      Endif  
      
201   Continue  
 
C---------------------------------------------------------------------------C
C             Set the normal vector of segment                              C
C---------------------------------------------------------------------------C      
 
      If (NSTEP.EQ.ISTART) Then
      Do N = 1,NUME2(NB)
          PDMIN = 10000
          NPMIN = 0
          BXN(NB,N,4) = 0.5 * (BXY2(N,1,1) + BXY2(N,2,1))
          BXN(NB,N,5) = 0.5 * (BXY2(N,1,2) + BXY2(N,2,2))
          BXN(NB,N,6) = 0
          Do I = 1,IJM
          If (IDC(I,1).EQ.1) Then
              PDIS = (CXY(I,1) - BXN(NB,N,4))**2 +
     &               (CXY(I,2) - BXN(NB,N,5))**2 
              If (PDIS.LE.PDMIN) Then
                  PDMIN = PDIS
                  NPMIN = I
              Endif
          Endif
          Enddo
          
          PA = BXY2(N,2,2) - BXY2(N,1,2)
          PB = BXY2(N,1,1) - BXY2(N,2,1)
          N_MAG = Sqrt(PA**2 + PB**2)
          PA = PA / N_MAG
          PB = PB / N_MAG
          N_SCA = PA * (CXY(NPMIN,1) - BXN(NB,N,4)) +
     &            PB * (CXY(NPMIN,2) - BXN(NB,N,5))
          If (N_SCA.LE.0) Then
              BXN(NB,N,1) = PA
              BXN(NB,N,2) = PB
              BXN(NB,N,3) = 0
          Else if (N_SCA.GT.0) Then
              BXN(NB,N,1) = -1 * PA
              BXN(NB,N,2) = -1 * PB
              BXN(NB,N,3) = 0
          Endif
          RN0(NB,N,1) = BXN(NB,N,1)
          RN0(NB,N,2) = BXN(NB,N,2)
          RN0(NB,N,3) = BXN(NB,N,3)
      Enddo
      Endif
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      
200   Continue
                               
      Endif
      
C===========================================================================c     
C         3D Immersed boundary detection     
C===========================================================================c

      If (IBMDOMAIN.EQ.3) Then
          
      Do 300 NB = 1,NUMIB
          
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C  
      Do I = 1,10000
      Do K = 1,3
          BXY3(I,K,1) = 0
          BXY3(I,K,2) = 0
          BXY3(I,K,3) = 0
      Enddo
      Enddo
      
      Do I = 1, NUME3(NB)
          BXY3(I,1,1) = BXY(NB,IBE3(NB,I,1),1)
          BXY3(I,1,2) = BXY(NB,IBE3(NB,I,1),2)
          BXY3(I,1,3) = BXY(NB,IBE3(NB,I,1),3)

          BXY3(I,2,1) = BXY(NB,IBE3(NB,I,2),1)
          BXY3(I,2,2) = BXY(NB,IBE3(NB,I,2),2)
          BXY3(I,2,3) = BXY(NB,IBE3(NB,I,2),3)

          BXY3(I,3,1) = BXY(NB,IBE3(NB,I,3),1)
          BXY3(I,3,2) = BXY(NB,IBE3(NB,I,3),2)
          BXY3(I,3,3) = BXY(NB,IBE3(NB,I,3),3)
      Enddo

! Get XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN
      XMAX = BXY(NB,1,1)
      XMIN = BXY(NB,1,1)
      YMAX = BXY(NB,1,2)
      YMIN = BXY(NB,1,2)
      ZMAX = BXY(NB,1,3)
      ZMIN = BXY(NB,1,3)

      Do I = 1,IBN(NB)
          If (BXY(NB,I,1).GE.XMAX) Then
              XMAX = BXY(NB,I,1)
          Endif
          If (BXY(NB,I,1).LE.XMIN) Then
              XMIN = BXY(NB,I,1)
          Endif
          If (BXY(NB,I,2).GE.YMAX) Then
              YMAX = BXY(NB,I,2)
          Endif
          If (BXY(NB,I,2).LE.YMIN) Then
              YMIN = BXY(NB,I,2)
          Endif
          If (BXY(NB,I,3).GE.ZMAX) Then
              ZMAX = BXY(NB,I,3)
          Endif
          If (BXY(NB,I,3).LE.ZMIN) Then
              ZMIN = BXY(NB,I,3)
          Endif
      Enddo 
      
! Solid cell detection 
      Do K = 1,KBM   
      Do 301 I = 1,IJM   
C---------------------------------------------------------------------------C
C             Ray-tracing Method                                            C
C---------------------------------------------------------------------------C             
      If ((CXY(I,1).LE.XMAX).AND.(CXY(I,1).GE.XMIN).AND. 
     &    (CXY(I,2).LE.YMAX).AND.(CXY(I,2).GE.YMIN).AND.
     &    (ZRC(I,K).LE.ZMAX).AND.(ZRC(I,K).GE.ZMIN)) Then 
      
      S1 = 0
      Do 302 N = 1,NUME3(NB)
      GK = (BXY3(N,2,1)-BXY3(N,1,1))*(BXY3(N,3,2)-BXY3(N,1,2))
     &    -(BXY3(N,2,2)-BXY3(N,1,2))*(BXY3(N,3,1)-BXY3(N,1,1))
      If (GK.NE.0) Then
          is_coincide = (((CXY(I,1).EQ.BXY3(N,1,1)).AND.
     &                (CXY(I,2).EQ.BXY3(N,1,2))).OR.
     &                ((CXY(I,1).EQ.BXY3(N,2,1)).AND.
     &                (CXY(I,2).EQ.BXY3(N,2,2))).OR.
     &                ((CXY(I,1).EQ.BXY3(N,3,1)).AND.
     &                (CXY(I,2).EQ.BXY3(N,3,2))))
          If (is_coincide) Then
              S1 = 999
              Exit
          Endif

          CRP(1) = (BXY3(N,1,1)-CXY(I,1))*(BXY3(N,2,2)-CXY(I,2))
     &            -(BXY3(N,1,2)-CXY(I,2))*(BXY3(N,2,1)-CXY(I,1))
          CRP(2) = (BXY3(N,2,1)-CXY(I,1))*(BXY3(N,3,2)-CXY(I,2))
     &            -(BXY3(N,2,2)-CXY(I,2))*(BXY3(N,3,1)-CXY(I,1))
          CRP(3) = (BXY3(N,3,1)-CXY(I,1))*(BXY3(N,1,2)-CXY(I,2))
     &            -(BXY3(N,3,2)-CXY(I,2))*(BXY3(N,1,1)-CXY(I,1))

          If (((CRP(1).GE.0).AND.(CRP(2).GE.0).AND.(CRP(3).GE.0)).OR.
     &        ((CRP(1).LT.0).AND.(CRP(2).LT.0).AND.(CRP(3).LT.0))) Then

          FL=((BXY3(N,2,1)-BXY3(N,1,1))*(BXY3(N,3,2)-BXY3(N,1,2))-
     &        (BXY3(N,2,2)-BXY3(N,1,2))*(BXY3(N,3,1)-BXY3(N,1,1))) 
          FR =(((CXY(I,2)-BXY3(N,1,2))*(BXY3(N,2,1)-BXY3(N,1,1))*
     &        (BXY3(N,3,3)-BXY3(N,1,3)))+
     &        ((CXY(I,1)-BXY3(N,1,1))*(BXY3(N,2,3)-BXY3(N,1,3))*
     &        (BXY3(N,3,2)-BXY3(N,1,2)))-
     &        ((CXY(I,1)-BXY3(N,1,1))*(BXY3(N,2,2)-BXY3(N,1,2))*
     &        (BXY3(N,3,3)-BXY3(N,1,3)))-
     &        ((CXY(I,2)-BXY3(N,1,2))*(BXY3(N,2,3)-BXY3(N,1,3))*
     &        (BXY3(N,3,1)-BXY3(N,1,1))))
          ZIB = FR/FL+BXY3(N,1,3)
          
          If (ZRC(I,K).LT.ZIB) Then
              S1 = S1 + 1
          Else If (ZRC(I,K).EQ.ZIB) Then
              S1 = 1
              Exit
          Endif
          Endif
      Endif 
                   
302   Continue
      
      If (mod(S1,2) .GT. 0.0) Then
          IDC(I,K) = 1
      Endif
      Endif               

301   Continue 
      Enddo 

C---------------------------------------------------------------------------C
C             Set the normal vector of segment                              C
C---------------------------------------------------------------------------C 
      
      If (NSTEP.EQ.ISTART) Then
      Do N = 1,NUME3(NB)
          PDMIN = 10000
          NPMIN = 0
          BXN(NB,N,4) = (BXY3(N,1,1)+BXY3(N,2,1)+BXY3(N,3,1)) / 3
          BXN(NB,N,5) = (BXY3(N,1,2)+BXY3(N,2,2)+BXY3(N,3,2)) / 3
          BXN(NB,N,6) = (BXY3(N,1,3)+BXY3(N,2,3)+BXY3(N,3,3)) / 3
          Do K = 1,KBM
          Do I = 1,IJM
          If (IDC(I,K).EQ.1) Then
              PDIS = (CXY(I,1) - BXN(NB,N,4))**2 +
     &               (CXY(I,2) - BXN(NB,N,5))**2 + 
     &               (ZRC(I,K) - BXN(NB,N,6))**2 
              If (PDIS.LE.PDMIN) Then
                  PDMIN = PDIS
                  NPMIN = I
              Endif
          Endif
          Enddo
          Enddo
          PA = (BXY3(N,2,2)-BXY3(N,1,2)) * (BXY3(N,3,3)-BXY3(N,1,3))
     &        -(BXY3(N,2,3)-BXY3(N,1,3)) * (BXY3(N,3,2)-BXY3(N,1,2))
          PB = (BXY3(N,2,3)-BXY3(N,1,3)) * (BXY3(N,3,1)-BXY3(N,1,1))
     &        -(BXY3(N,2,1)-BXY3(N,1,1)) * (BXY3(N,3,3)-BXY3(N,1,3))
          PC = (BXY3(N,2,1)-BXY3(N,1,1)) * (BXY3(N,3,2)-BXY3(N,1,2))
     &        -(BXY3(N,2,2)-BXY3(N,1,2)) * (BXY3(N,3,1)-BXY3(N,1,1))
          PD = -1 * (PA*BXY(N,1,1) + PB*BXY(N,1,2) + PC*BXY(N,1,3))
          
          N_MAG = Sqrt(PA**2 + PB**2 + PC**2)
          PA = PA / N_MAG
          PB = PB / N_MAG
          PC = PC / N_MAG
          N_SCA = PA * (CXY(NPMIN,1) - BXN(NB,N,4)) +
     &            PB * (CXY(NPMIN,2) - BXN(NB,N,5)) +
     &            PC * (ZRC(NPMIN,K) - BXN(NB,N,6))
          If (N_SCA.LE.0) Then
              BXN(NB,N,1) = PA
              BXN(NB,N,2) = PB
              BXN(NB,N,3) = PC
          Else if (N_SCA.GT.0) Then
              BXN(NB,N,1) = -1 * PA
              BXN(NB,N,2) = -1 * PB
              BXN(NB,N,3) = -1 * PC
          Endif
          RN0(NB,N,1) = BXN(NB,N,1)
          RN0(NB,N,2) = BXN(NB,N,2)
          RN0(NB,N,3) = BXN(NB,N,3)
      Enddo
      Endif
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C  

300   Continue

      Endif           
      
C===========================================================================C
C         Classify ghost cells and solid cells
C===========================================================================C 
          
      Do K = 1,KBM
      Do I = 1,IJM
          If (IDC(I,K).LT.0.5) Then
              S1 = 0
              Do J = 1,CELL_POLYGEN(I)
                  If (IDC(CELL_SIDE(I,J,2),K).EQ.1) Then
                      S1 = S1 + 1
                  Endif
              Enddo
              If (S1.EQ.CELL_POLYGEN(I)) IDC(I,K) = 1   
          Endif
      Enddo
      Enddo
      
      Do K = 1,KBM
      Do I = 1,IJM
          If (IDC(I,K).EQ.1) Then
              S1 = 0
              Do J = 1,CELL_POLYGEN(I)
                  If (IDC(CELL_SIDE(I,J,2),K).LT.0.5) Then
                      S1 = S1 + 1
                  Endif
              Enddo
              If (S1.EQ.CELL_POLYGEN(I)) IDC(I,K) = 0 
          Endif
      Enddo
      Enddo
      
      Do K = 1,KBM
      Do I = 1,IJM
          If (IDC(I,K).EQ.1) Then  
              Do J = 1,CELL_POLYGEN(I)
              If ((CELL_SIDE(I,J,2).GT.0).AND.
     &           (IDC(CELL_SIDE(I,J,2),K).LT.0.9)) Then
                  IDC(CELL_SIDE(I,J,2),K) = 2
              Endif
              If (K.EQ.1) Then
                  If (IDC(I,2).NE.1) IDC(I,2) = 3
              Else if (K.EQ.KBM) Then
                  If (IDC(I,KBM-1).NE.1) IDC(I,KBM-1) = 4
              Else
                  If (IDC(I,K-1).NE.1) IDC(I,K-1) = 4
                  If (IDC(I,K+1).NE.1) IDC(I,K+1) = 3
              Endif
              Enddo
          
              Do J = 1,CELL_POLYGEN(I)
              Do N = 1,POINT_CELL(CELL_SIDE(I,J,3))
              If ((INDEX_POINT(CELL_SIDE(I,J,3),N).GT.0).AND.
     &           (IDC(INDEX_POINT(CELL_SIDE(I,J,3),N),K).LT.0.1)) Then
                  IDC(INDEX_POINT(CELL_SIDE(I,J,3),N),K) = 2
              Endif
              Enddo
              Enddo
          Endif
      Enddo
      Enddo
      
      Do K = 1,KBM
      Do I = 1,IJM
          If (IDC(I,K).EQ.4) Then 
              Do K1 = K+1, KBM
                  IDE(I,K1) = 1
              Enddo
          Endif
          If (IDC(I,K).EQ.1) IDE(I,K) = 1
          If ((IFC(I,K).GE.1).AND.(IDC(I,K).NE.1)) Then
              IFC(I,K) = 1
          Else
              IFC(I,K) = 0
          Endif
      Enddo
      Enddo
      
      

C===========================================================================C
C         IBM cells variable assignment
C===========================================================================C 

      N_SL = 0
      N_GH = 0
     
      Do K = 1,KBM
      Do I = 1,IJM
          If (IDC(I,K).EQ.1) Then
              N_SL = N_SL + 1
              ISC(N_SL,1) = I
              ISC(N_SL,2) = K
          Endif
          If (IDC(I,K).GT.1.9) Then
              N_GH = N_GH + 1
              IGC(N_GH,1) = I
              IGC(N_GH,2) = K 
          Endif
      Enddo
      Enddo
      
      NUMSL = N_SL
      NUMGH = N_GH
      
      ! search image point
      If (IBMDOMAIN.EQ.2) Then
      Do N = 1,NUMGH
          
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      PDMIN = 10000
      NPMIN = 1
      NBMIN = 1
      
      Do I = 1,2
      Do J = 1,2
          LINE(I,J) = 1
      Enddo
      Enddo
	               
      Do NB = 1,NUMIB
      Do I = 1, NUME2(NB)
          P1 = IBE2(NB,I,1)
          P2 = IBE2(NB,I,2)
          PDIS = (CXY(IGC(N,1),1) - 0.5*(BXY(NB,P1,1)+BXY(NB,P2,1)))**2
     &          +(CXY(IGC(N,1),2) - 0.5*(BXY(NB,P1,2)+BXY(NB,P2,2)))**2
          
          If (PDIS.LE.PDMIN) Then 
              PDMIN = PDIS
              NPMIN = I
              NBMIN = NB
              LINE(1,1) = BXY(NB,P1,1)
              LINE(1,2) = BXY(NB,P1,2)
              LINE(2,1) = BXY(NB,P2,1)
              LINE(2,2) = BXY(NB,P2,2)
          Endif
      Enddo
      Enddo
      
      If (IDC(IGC(N,1),IGC(N,2)).EQ.2) Then
          PA = LINE(2,2)-LINE(1,2)
          PB = LINE(1,1)-LINE(2,1)
          PC = LINE(1,2)*LINE(2,1)-LINE(2,2)*LINE(1,1)
          XIM = (-2*PA*PC+(PB**2-PA**2)*CXY(IGC(N,1),1)
     &           -2*PA*PB*CXY(IGC(N,1),2)) / (PA**2+PB**2)
          YIM = (-2*PB*PC+(PA**2-PB**2)*CXY(IGC(N,1),2)
     &           -2*PA*PB*CXY(IGC(N,1),1)) / (PA**2+PB**2)
          ZIM = ZRC(IGC(N,1),IGC(N,2))
          
          IPX = 1.5*CXY(IGC(N,1),1) - 0.5*XIM
          IPY = 1.5*CXY(IGC(N,1),2) - 0.5*YIM
          IPZ = ZIM
          IPK = IGC(N,2)
          IPNV(1) = BXN(NBMIN,NPMIN,1)
          IPNV(2) = BXN(NBMIN,NPMIN,2)
          IPNV(3) = 0
          IPDIS = Abs(PA*IPX+PB*IPY+PC) / Sqrt(PA**2+PB**2)
      Else if (IDC(IGC(N,1),IGC(N,2)).GE.3) Then
          XIM = CXY(IGC(N,1),1)
          YIM = CXY(IGC(N,1),2)
          IPX = XIM
          IPY = YIM
          If (IDC(IGC(N,1),IGC(N,2)).EQ.3) Then
              If (IGC(N,2).EQ.1) IPK = 1
              If (IGC(N,2).GT.1) IPK = IGC(N,2) + 1
              IPZ = ZRC(IGC(N,1),IPK)
              IPNV(3) = -1
              IPDIS = Abs(IPZ - Z2D(1))   
          Else if (IDC(IGC(N,1),IGC(N,2)).EQ.4) Then
              If (IGC(N,2).EQ.KBM) IPK = KBM
              If (IGC(N,2).LT.KBM) IPK = IGC(N,2) - 1
              IPZ = ZRC(IGC(N,1),IPK)
              IPNV(3) = 1
              IPDIS = Abs(IPZ - Z2D(2))
          Endif  
          IPZ = ZRC(IGC(N,1),IPK)
          IPNV(1) = 0
          IPNV(2) = 0
      Endif
     
      BGP(N,1) = IPX
      BGP(N,2) = IPY
      BGP(N,3) = IPZ
      BGP(N,4) = IPK
      BGP(N,5) = NBMIN
      BGP(N,6) = NPMIN
      BGP(N,7) = IPNV(1)
      BGP(N,8) = IPNV(2)
      BGP(N,9) = IPNV(3)
      BGP(N,10) = IPDIS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
          
      Enddo
      Endif
      
      If (IBMDOMAIN.EQ.3) Then
      Do N = 1,NUMGH
          
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      PDMIN = 10000
      NPMIN = 1
      NBMIN = 1
      
      Do I = 1,3
      Do J = 1,3
          PLANE(I,J) = 1
      Enddo
      Enddo
            
      Do NB = 1,NUMIB
      Do I = 1, NUME3(NB)
          P1 = IBE3(NB,I,1)
          P2 = IBE3(NB,I,2)
          P3 = IBE3(NB,I,3)
          PDIS = (CXY(IGC(N,1),1) - (BXY(NB,P1,1) +
     &            BXY(NB,P2,1) + BXY(NB,P3,1)) / 3)**2 + 
     &           (CXY(IGC(N,1),2) - (BXY(NB,P1,2) +
     &            BXY(NB,P2,2) + BXY(NB,P3,2)) / 3)**2 +
     &           (ZRC(IGC(N,1),IGC(N,2)) - (BXY(NB,P1,3) +
     &            BXY(NB,P2,3) + BXY(NB,P3,3)) / 3)**2 

          If (PDIS.LE.PDMIN) Then 
              PDMIN = PDIS
              NPMIN = I
              NBMIN = NB
              PLANE(1,1) = BXY(NB,P1,1)
              PLANE(1,2) = BXY(NB,P1,2)
              PLANE(1,3) = BXY(NB,P1,3)
              PLANE(2,1) = BXY(NB,P2,1)
              PLANE(2,2) = BXY(NB,P2,2)
              PLANE(2,3) = BXY(NB,P2,3)
              PLANE(3,1) = BXY(NB,P3,1)
              PLANE(3,2) = BXY(NB,P3,1)
              PLANE(3,3) = BXY(NB,P3,1)
          Endif
      Enddo
      Enddo
      
      PA = (PLANE(2,2)-PLANE(1,2)) * (PLANE(3,3)-PLANE(1,3))
     &    -(PLANE(2,3)-PLANE(1,3)) * (PLANE(3,2)-PLANE(1,2))
      PB = (PLANE(2,3)-PLANE(1,3)) * (PLANE(3,1)-PLANE(1,1))
     &    -(PLANE(2,1)-PLANE(1,1)) * (PLANE(3,3)-PLANE(1,3))
      PC = (PLANE(2,1)-PLANE(1,1)) * (PLANE(3,2)-PLANE(1,2))
     &    -(PLANE(2,2)-PLANE(1,2)) * (PLANE(3,1)-PLANE(1,1))
      PD = -1 * (PA*PLANE(1,1)+PB*PLANE(1,2)+PC*PLANE(1,3))
      
      TLINE = -1 * (PA*CXY(IGC(N,1),1) + PB*CXY(IGC(N,1),2)+
     &         PC*ZRC(IGC(N,1),IGC(N,2)) + PD) /
     &        (PA**2 + PB**2 + PC**2)
     
      XIM = 2*PA*TLINE + CXY(IGC(N,1),1)
      YIM = 2*PB*TLINE + CXY(IGC(N,1),2)
      ZIM = 2*PC*TLINE + ZRC(IGC(N,1),IGC(N,2))
      
      IPX = 1.5*CXY(IGC(N,1),1) - 0.5*XIM
      IPY = 1.5*CXY(IGC(N,1),2) - 0.5*YIM
      IPZ = 1.5*ZRC(IGC(N,1),IGC(N,2)) - 0.5*ZIM
      IPNV(1) = BXN(NBMIN,NPMIN,1)
      IPNV(2) = BXN(NBMIN,NPMIN,2)
      IPNV(3) = BXN(NBMIN,NPMIN,3)
      IPDIS = Abs(PA*XIM+PB*YIM+PC*ZIM+PD) / Sqrt(PA**2+PB**2+PC**2)
      
      ZNE = 1
      Do While (ZRC(IGC(N,1),ZNE).GE.IPZ)
          ZNE = ZNE + 1
      Enddo
      IPK = ZNE

      BGP(N,1) = IPX
      BGP(N,2) = IPY
      BGP(N,3) = IPZ
      BGP(N,4) = IPK
      BGP(N,5) = NBMIN
      BGP(N,6) = NPMIN
      BGP(N,7) = IPNV(1)
      BGP(N,8) = IPNV(2)
      BGP(N,9) = IPNV(3)
      BGP(N,10) = IPDIS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      
      Enddo
      Endif
      
      End Subroutine
