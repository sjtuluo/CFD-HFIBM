C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (21/02/2023)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                                                                           #     
C                     Visualization of IBM results                          #                                           
C                                                                           #
C############################################################################
      Subroutine IBMVISUAL(FIGMOD,FRAMENUM)

      Include './Include/OCERM_INF'
      
      Integer FIGMOD,FRAMENUM
      
      
      If (IBMSCHEME.GE.0) Then
C===========================================================================C     
C         Visualization of ghost cells  
C===========================================================================C  
      
      If (FIGMOD.EQ.1) Then
          Do NB = 1,NUMIB
          Do J = 1,IBN(NB)   
              Write(IBMB,1000) BXY(NB,J,1),BXY(NB,J,2),BXY(NB,J,3)
          Enddo
          If (IBMDOMAIN.EQ.2) Then 
              Do I = 1,NUME2(NB)
                  Write(IBMN,'(4F16.8)') BXN(NB,I,4),BXN(NB,I,5),
     &                                   BXN(NB,I,1),BXN(NB,I,2)
              Enddo
          Else if (IBMDOMAIN.EQ.3) Then 
              Do I = 1,NUME3(NB)
                  Write(IBMN,'(6F16.8)') BXN(NB,I,4),BXN(NB,I,5),
     &            BXN(NB,I,6),BXN(NB,I,1),BXN(NB,I,2),BXN(NB,I,3)
                  Write(IBMOBJ,'(3I8,3F16.8)')IBE3(NB,I,1),IBE3(NB,I,2),
     &            IBE3(NB,I,3),BXN(NB,I,1),BXN(NB,I,2),BXN(NB,I,3)
              Enddo
          Endif
          Enddo

          
          Do K = 1,KBM
          Do I = 1,IJM
              If ((IDC(I,K).LT.0.1)) Then
                  Write(IBMF,1000) CXY(I,1),CXY(I,2),ZRC(I,K)
              Else if ((IDC(I,K).EQ.1)) Then
                  Write(IBMS,1000) CXY(I,1),CXY(I,2),ZRC(I,K)
              Else If (IDC(I,K).GT.1.9) Then
                  Write(IBMGC,1000) CXY(I,1),CXY(I,2),ZRC(I,K)
              Endif
          Enddo
          Enddo
          
          If (IBGCSCHEME.GT.0) Then
          near_points = 28
          Do N = 1,NUMGH
              Write(IBMGCI,1000) BGP(N,1),BGP(N,2),BGP(N,3)
              Do J = 1,near_points
              I = NBR(IGC(N,1),J)
              K = BGP(N,4)
              If (IDC(I,K).LT.0.9) Then
                  DK = (BGP(N,1)-CXY(I,1))**2 + (BGP(N,2)-CXY(I,2))**2
                  If (DK.LE.RITP**2) Then
                  Write(IBMP,1000) CXY(I,1),CXY(I,2),ZRC(I,K)
                  Endif
              Endif
              Enddo
          Enddo
          Endif
          
      Endif
      
C===========================================================================C      
C         Visualization of ghost points      
C===========================================================================C  

      If (FIGMOD.EQ.2) Then        
          Do NB = 1,NUMIB
          Do J = 1,IBN(NB)   
              Write(IBMB,1000) BXY(NB,J,1),BXY(NB,J,2),BXY(NB,J,3)
          Enddo
          Enddo
      
          Do K = 1,KB
          Do I = 1,IJP
              If (IDP(I,K).LT.0.5) Then
                  Write(IBMF,1000) PXY(I,1),PXY(I,2),ZRP(I,K)
              Else if (IDP(I,K).GT.0.5) Then
                  Write(IBMS,1000) PXY(I,1),PXY(I,2),ZRP(I,K)
              Endif
          Enddo
          Enddo
      Endif

C===========================================================================C
C         Visualization of IBM mesh
C===========================================================================C

      If (FIGMOD.EQ.3) Then
          Do NB = 1,NUMIB
          Do J = 1,IBN(NB)   
              Write(IBMB,1000) BXY(NB,J,1),BXY(NB,J,2),BXY(NB,J,3)
          Enddo
          Enddo

          Do I = 1,IJM
              If((CXY(I,1).GT.-0.1).AND.(CXY(I,1).LT.0.1).AND.
     &           (CXY(I,2).GT.-0.1).AND.(CXY(I,2).LT.0.1)) Then
          Do J = 1,CELL_POLYGEN(I)
              Write(IBMF,1000) PXY(CELL_SIDE(I,J,3),1),
     &                         PXY(CELL_SIDE(I,J,3),2),ZRC(I,1)
          Enddo
          Write(IBMF,1000) PXY(CELL_SIDE(I,1,3),1),
     &                     PXY(CELL_SIDE(I,1,3),2),ZRC(I,1)
          Endif
          Enddo
          
          K = 5
          Do I = 1,IJM
              If ((IDC(I,K).EQ.1)) Then
              Do J = 1,CELL_POLYGEN(I)
                  Write(IBMS,1000) PXY(CELL_SIDE(I,J,3),1),
     &                             PXY(CELL_SIDE(I,J,3),2),ZRC(I,K)
              Enddo
              Else If (IDC(I,K).GT.1.9) Then
              Do J = 1,CELL_POLYGEN(I)
                  Write(IBMGC,1000) PXY(CELL_SIDE(I,J,3),1),
     &                              PXY(CELL_SIDE(I,J,3),2),ZRC(I,K)
              Enddo
              Endif
          Enddo
      Endif

C===========================================================================C
C         Visualization of IBM moving boundary
C===========================================================================C

      If (FIGMOD.EQ.4) Then
           If (NSTEP.EQ.100) Then
              Do I = 1,NUME3(1)
              BXN(2,I,1) = BXN(1,I,1)
              BXN(2,I,2) = BXN(1,I,2)
              BXN(2,I,3) = BXN(1,I,3)
              Enddo
           Endif
           
          If (Mod(NSTEP,FRAMENUM) .EQ. 0.0) Then
              
          Write(IBMB,1100) 0,0,100
          Do I = 1,IBN(1)
              Write(IBMB,1000) BXY(1,I,1),BXY(1,I,2),BXY(1,I,3)
          Enddo
          
          Write(IBMDIB,1100) 0,0,100
          Do J = 1,IBN(2)   
              Write(IBMDIB,1000) BXY(2,J,1),BXY(2,J,2),BXY(2,J,3)
          Enddo
      
          Write(IBMS,1100) 0,0,100
          Do K = 1,KBM
          Do I = 1,IJM
              If(IDC(I,K).EQ.1) THEN
                  Write(IBMS,1000) CXY(I,1),CXY(I,2),ZRC(I,K)
              Endif
          Enddo
          Enddo
          
          Write(IBMN,1100) 0,0,100
          Write(IBMT,1100) 0,0,100
          NB = 1
          If (IBMDOMAIN.EQ.2) Then 
              Do I = 1,NUME2(NB)
                  Write(IBMN,'(2F16.8)') BXN(NB,I,1),BXN(NB,I,2) 
              Enddo
          Else if (IBMDOMAIN.EQ.3) Then 
              Do I = 1,NUME3(NB)
                  Write(IBMN,1000) BXN(NB,I,1),BXN(NB,I,2),BXN(NB,I,3)
                  Write(IBMT,1000) BXN(2,I,1),BXN(2,I,2),BXN(2,I,3)
              Enddo
          Endif
          Endif
      Endif
          
      Endif

1000  Format (3F16.8) 
1100  Format (3I16)
      
      
      End Subroutine