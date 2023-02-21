C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (21/02/2023)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                                                                           #
C                Detecte differnet kind of points.                          #                                           
C                                                                           #
C############################################################################      
      Subroutine IBMALIDP

      Include './Include/OCERM_INF'
      
      Real(kind=8) BXY2(10000,2,2),BXY3(10000,3,3),CRP(3),Z2D(2)
      Double precision XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN,FL,FR,ZIB,GK
      Integer S1
      

C===========================================================================C      
C         Initialize the arrays      
C===========================================================================C

      Z2D(1) = -999
      Z2D(2) = -0.05
      
      Do K = 1,KB
      Do I = 1,IJP
          IDP(I,K) = 0
          ZRP(I,K) = Z(K)*(HP(I)+ELFV(I))+ELFV(I)
      Enddo
      Enddo

      Do I = 1,3
          CRP(I) = 0
      Enddo
      
C===========================================================================C      
C      2D Immersed boundary detection     
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
      Do 201 I = 1,IJP              
C---------------------------------------------------------------------------C
C             Ray-tracing Method                                            C
C---------------------------------------------------------------------------C           
      If ((PXY(I,1).LT.XMAX).AND.(PXY(I,1).GT.XMIN).AND. 
     &    (PXY(I,2).LT.YMAX).AND.(PXY(I,2).GT.YMIN)) Then
     &
      S1 = 0   
      
      Do 202 N = 1,NUME2(NB)
          CRP(1) = (PXY(I,1)-BXY2(N,1,1)) * (BXY2(N,2,2)-PXY(I,2))
     &            -(BXY2(N,2,1)-PXY(I,1)) * (PXY(I,2)-BXY2(N,1,2))
                           
          If (((BXY2(N,1,1).LE.PXY(I,1)).AND.
     &        (BXY2(N,2,1).GT.PXY(I,1)).AND.(CRP(1).GE.0))
     &        .OR.
     &        ((BXY2(N,1,1).GE.PXY(I,1)).AND.
     &        (BXY2(N,2,1).LT.PXY(I,1)).AND.(CRP(1).LE.0))) 
     &    Then               
              S1 = S1 + 1
              If (CRP(1) .EQ. 0) Then     
                  S1 = 1
                  Exit
              Endif
          Endif  
202   Continue  
      
      Do K = 1,KB
      If ((ZRP(I,K).GE.Z2D(1)).AND.(ZRP(I,K).LE.Z2D(2))) Then
          If (mod(S1,2) .GT. 0) Then
              IDP(I,K) = 1     
          Endif
      Endif
      Enddo
      
      Endif

201   Continue
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      
200   Continue

      Endif 
      
C===========================================================================C      
C        3D Immersed boundary detection      
C===========================================================================C

      If (IBMDOMAIN.EQ.3) Then
      
      Do 300 NB = 1,NUMIB
          
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C   
      Do I = 1,10000
      Do K = 1,3
          BXY3(I,K,1) = 0
          BXY3(I,K,1) = 0
          BXY3(I,K,1) = 0
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
      Do K = 1,KB    
      Do 301 I = 1,IJP    

C***************************************************************************C
C             Ray-casting Method Begins                                     C
C***************************************************************************C            
      If ((PXY(I,1).LE.XMAX).AND.(PXY(I,1).GE.XMIN).AND. 
     &    (PXY(I,2).LE.YMAX).AND.(PXY(I,2).GE.YMIN).AND.
     &    (ZRP(I,K).LE.ZMAX).AND.(ZRP(I,K).GE.ZMIN)) Then 
     
      S1 = 0
      Do 302 N = 1,NUME3(NB)
      GK = (BXY3(N,2,1)-BXY3(N,1,1))*(BXY3(N,3,2)-BXY3(N,1,2))
     &    -(BXY3(N,2,2)-BXY3(N,1,2))*(BXY3(N,3,1)-BXY3(N,1,1))
      If (GK.NE.0) Then
          is_coincide = (((PXY(I,1).EQ.BXY3(N,1,1)).AND.
     &                (PXY(I,2).EQ.BXY3(N,1,2))).OR.
     &                ((PXY(I,1).EQ.BXY3(N,2,1)).AND.
     &                (PXY(I,2).EQ.BXY3(N,2,2))).OR.
     &                ((PXY(I,1).EQ.BXY3(N,3,1)).AND.
     &                (PXY(I,2).EQ.BXY3(N,3,2))))
          If (is_coincide) Exit

          CRP(1) = (BXY3(N,1,1)-PXY(I,1))*(BXY3(N,2,2)-PXY(I,2))
     &            -(BXY3(N,1,2)-PXY(I,2))*(BXY3(N,2,1)-PXY(I,1))
          CRP(2) = (BXY3(N,2,1)-PXY(I,1))*(BXY3(N,3,2)-PXY(I,2))
     &            -(BXY3(N,2,2)-PXY(I,2))*(BXY3(N,3,1)-PXY(I,1))
          CRP(3) = (BXY3(N,3,1)-PXY(I,1))*(BXY3(N,1,2)-PXY(I,2))
     &            -(BXY3(N,3,2)-PXY(I,2))*(BXY3(N,1,1)-PXY(I,1))

          If (((CRP(1).GE.0).AND.(CRP(2).GE.0).AND.(CRP(3).GE.0)).OR.
     &        ((CRP(1).LE.0).AND.(CRP(2).LE.0).AND.(CRP(3).LE.0))) Then
          If ((CRP(1).EQ.0).OR.(CRP(2).EQ.0).OR.(CRP(3).EQ.0)) Exit

          FL=((BXY3(N,2,1)-BXY3(N,1,1))*(BXY3(N,3,2)-BXY3(N,1,2))-
     &        (BXY3(N,2,2)-BXY3(N,1,2))*(BXY3(N,3,1)-BXY3(N,1,1))) 
          FR =(((PXY(I,2)-BXY3(N,1,2))*(BXY3(N,2,1)-BXY3(N,1,1))*
     &        (BXY3(N,3,3)-BXY3(N,1,3)))+
     &        ((PXY(I,1)-BXY3(N,1,1))*(BXY3(N,2,3)-BXY3(N,1,3))*
     &        (BXY3(N,3,2)-BXY3(N,1,2)))-
     &        ((PXY(I,1)-BXY3(N,1,1))*(BXY3(N,2,2)-BXY3(N,1,2))*
     &        (BXY3(N,3,3)-BXY3(N,1,3)))-
     &        ((PXY(I,2)-BXY3(N,1,2))*(BXY3(N,2,3)-BXY3(N,1,3))*
     &        (BXY3(N,3,1)-BXY3(N,1,1))))
          ZIB = FR/FL+BXY3(N,1,3)
          
          If (ZRP(I,K).LT.ZIB) Then
              S1 = S1 + 1
          Else If (ZRP(I,K).EQ.ZIB) Then
              S1 = 1
              Exit
          Endif
          Endif
      Endif 
                   
302   Continue
      
      If (mod(S1,2) .GT. 0) Then
          IDP(I,K) = 1
      Endif
      Endif
C***************************************************************************C
C             End of Ray-casting method                                     C
C***************************************************************************C                  

301   Continue 
      
      Enddo
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      
300   Continue

      Endif  
      
      Do K = 1,KB
      Do I = 1,IJP
          If (IDP(I,K).LT.0.5) Then
              S1 = 0
              Do J = 1,POINT_CELL(I)
              If (K.EQ.KB) Then
                  If (IDC(INDEX_POINT(I,J),KBM).EQ.1) Then
                      S1 = S1 + 1
                  Endif
              Else
                  If (IDC(INDEX_POINT(I,J),K).EQ.1) Then
                      S1 = S1 + 1
                  Endif
              Endif
              Enddo
              If (S1.EQ.POINT_CELL(I)) IDP(I,K) = 1   
          Endif
      Enddo
      Enddo
      
      Do K = 1,KB
      Do I = 1,IJP
          If (IDP(I,K).EQ.1) Then
              S1 = 0
              Do J = 1,POINT_CELL(I)
              If (K.EQ.KB) Then
                  If (IDC(INDEX_POINT(I,J),KBM).LT.0.5) Then
                      S1 = S1 + 1
                  Endif
              Else
                  If (IDC(INDEX_POINT(I,J),K).LT.0.5) Then
                      S1 = S1 + 1
                  Endif
              Endif
              Enddo
              If (S1.EQ.POINT_CELL(I)) IDP(I,K) = 1   
          Endif
      Enddo
      Enddo

      
      End Subroutine