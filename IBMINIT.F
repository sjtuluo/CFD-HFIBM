C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (21/02/2023)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                                                                           #     
C      Read data of IBM and modify parameters of mesh geometry in 3D        #                                           
C                                                                           #
C############################################################################
      Subroutine IBMINIT
      
      Include './Include/OCERM_INF'
      
      Dimension COM(80)
      Real(Kind=8) X_SUM,Y_SUM,Z_SUM
      Integer ISNBR(IJM)
      Integer I2,I3,NUM
      Character*10 BOUNDARY_ID
      

C===========================================================================C      
C         Set the parameter of immersed boundary method 
C===========================================================================C

      Open (IBMSET,file='./Include/IBM/IBM_SETTING.dat')
      
      Do N = 1,3
          Read (IBMSET,'(100A1)') (COM(I),I = 1,80)
      Enddo
      
      Read (IBMSET,1000) (COM(I),I = 1,80)
      Write(IUPRT,1000) (COM(I),I = 1,80)
      Read (IBMSET,'(3I16)') IBMSCHEME,IBMDOMAIN,NUMIB
      Write(IUPRT,'(3I16)') IBMSCHEME,IBMDOMAIN,NUMIB
      Read (IBMSET,1000) (COM(I),I = 1,80)
      Write(IUPRT,1000) (COM(I),I = 1,80)
      Read (IBMSET,'(3I16)') IBFSCHEME,IBGCSCHEME,IBGCMODE
      Write(IUPRT,'(3I16)') IBFSCHEME,IBGCSCHEME,IBGCMODE
      Read (IBMSET,1000) (COM(I),I = 1,80)
      Write(IUPRT,1000) (COM(I),I = 1,80)
      Read (IBMSET,'(2F16.8,I16)') RITP,RIFL,INFCELL
      Write (IUPRT,'(2F16.8,I16)') RITP,RIFL,INFCELL
      
      Close (IBMSET)

C===========================================================================C      
C         Initialize the arrays 
C===========================================================================C

      Do K = 1,KBM
          Do J = 1,IJM
              IDC(I,K) = 0
              IFC(I,K) = 0
              IDE(I,K) = 0
              ZRC(I,K) = 0
          Enddo
      Enddo
      
      Do I = 1,IJM
          ISNBR(I) = 0
          Do N = 1,30
              NBR(I,N) = 0
          Enddo
      Enddo
      
      Do K = 1,KB
          Do I = 1,IJP
              IDP(I,K) = 0
              ZRP(I,K) = 0
		    SAU(I,K) = 0
	 	    SAV(I,K) = 0
		    SAW(I,K) = 0
          Enddo
      Enddo    
      
      Do NB = 1,5
          Do I = 1,IJM
              Do J = 1,3
                  BXY(NB,I,J) = 0
                  BXU(NB,I,J) = 0
                  RI0(NB,I,J) = 0
                  RN0(NB,I,J) = 0
              Enddo
              Do J = 1,2
                  IBE2(NB,I,J) = 0
                  IBE3(NB,I,J) = 0
                  BXN(NB,I,J) = 0
              Enddo
              IBE3(NB,I,3) = 0
              BXN(NB,I,3) = 1
          Enddo
          
          Do I = 1,3
              Do J = 1,3
                  RGIR(NB,I,J) = 0
              Enddo
          Enddo
          
          Do I = 1,3
              RGF(NB,3) = 0
              
              RGX(NB,I) = 0
              RGU(NB,I) = 0
              RGW(NB,I) = 0
              RGQT(NB,I) = 0
          Enddo
          
          RGQT(NB,4) = 1
          RGM(NB) = 0
          FRN(NB) = 0
          NUME2(NB) = 0
          NUME3(NB) = 0
      Enddo
      
      IBSTEP = 1
      
C===========================================================================C      
C         Get the boundary information of IBM 
C===========================================================================C

      If (IBMDOMAIN.EQ.2) Then 
          Open(IBM2D,File = './Include/IBM/IBM2D.dat')
          Do NB = 1,NUMIB     
              Read(IBM2D,*) BOUNDARY_ID       ! INDEX of boundary
              Read(IBM2D,*) IBN(NB)      ! Number of IBN points
              Do I = 1,IBN(NB)
                  Read(IBM2D,*) BXY(NB,I,1),BXY(NB,I,2)
              Enddo
              
              NUME2(NB) = IBN(NB)     ! NUME2 is the number of segment
              Do I = 1, IBN(NB)-1
                  IBE2(NB,I,1) = I
                  IBE2(NB,I,2) = I + 1
              Enddo
              IBE2(NB,IBN(NB),1) = IBN(NB) 
              IBE2(NB,IBN(NB),2) = 1
          Enddo
          Close(IBM2D)
          
          Do NB = 1,NUMIB
              X_SUM = 0
              Y_SUM = 0
              Do I = 1,IBN(NB)
                  X_SUM = X_SUM + BXY(NB,I,1)
                  Y_SUM = Y_SUM + BXY(NB,I,2)
              Enddo
              RGX(NB,1) = X_SUM / IBN(NB)
              RGX(NB,2) = Y_SUM / IBN(NB)
              Do I = 1,IBN(NB)
                  RI0(NB,I,1) = BXY(NB,I,1) - RGX(NB,1)
                  RI0(NB,I,2) = BXY(NB,I,2) - RGX(NB,2)
              Enddo
          Enddo
      Endif
          
      If (IBMDOMAIN.EQ.3) Then        
          Open(IBM3D,File = './Include/IBM/IBM3D.dat')
          Do NB = 1,NUMIB     ! Iterate each boundary
              Read(IBM3D,*) BOUNDARY_ID       ! Index of boundary
              Read(IBM3D,*) IBN(NB)      ! Number of IBN points
              Do I = 1,IBN(NB)
                  Read(IBM3D,*) BXY(NB,I,1),BXY(NB,I,2),BXY(NB,I,3)
              Enddo  
              
              NUME3(NB) = 0     ! NUME3 is the number of triangle element
              Do I = 1, IBN(NB)/3
                  NUME3(NB) = NUME3(NB) + 1
                  IBE3(NB,NUME3(NB),1) = 3*I - 2
                  IBE3(NB,NUME3(NB),2) = 3*I - 1
                  IBE3(NB,NUME3(NB),3) = 3*I
              Enddo
          Enddo
          Close(IBM3D)
          
          Do NB = 1,NUMIB
              X_SUM = 0
              Y_SUM = 0
              Z_SUM = 0
              Do I = 1,IBN(NB)
                  X_SUM = X_SUM + BXY(NB,I,1)
                  Y_SUM = Y_SUM + BXY(NB,I,2)
                  Z_SUM = Z_SUM + BXY(NB,I,3)
              Enddo
              RGX(NB,1) = X_SUM / IBN(NB)
              RGX(NB,2) = Y_SUM / IBN(NB)
              RGX(NB,3) = Z_SUM / IBN(NB)
              Do I = 1,IBN(NB)
                  RI0(NB,I,1) = BXY(NB,I,1) - RGX(NB,1)
                  RI0(NB,I,2) = BXY(NB,I,2) - RGX(NB,2)
                  RI0(NB,I,3) = BXY(NB,I,3) - RGX(NB,3)
              Enddo
          Enddo
      Endif
      
      Do I = 1,IBN(1)
          BXY(2,I,1) = BXY(1,I,1)
          BXY(2,I,2) = BXY(1,I,2)
          BXY(2,I,3) = BXY(1,I,3) - 0.4
          Write(IBMCP,'(3F16.8)') BXY(2,I,1),BXY(2,I,2),BXY(2,I,3)
      Enddo 

C===========================================================================C      
C         Set neighbor cells
C===========================================================================C
      
      Do I = 1,IJM
          ISNBR(I) = 1
	    NUM = 0
	    Do J = 1,CELL_POLYGEN(I)
              If(CELL_SIDE(I,J,2).GT.0) Then
                  ISNBR(CELL_SIDE(I,J,2)) = 1
		        NUM = NUM + 1
                  NBR(I,NUM) = CELL_SIDE(I,J,2)
              Endif
          Enddo
          
	    Do I2 = 1,IJM
		    If(ISNBR(I2).EQ.1) Then
		    Do J = 1,CELL_POLYGEN(I2)
			    If((CELL_SIDE(I2,J,2).GT.0).AND.
     &               (ISNBR(CELL_SIDE(I2,J,2)).EQ.0)) Then
				    ISNBR(CELL_SIDE(I2,J,2)) = 2
				    NUM = NUM + 1
                      NBR(I,NUM) = CELL_SIDE(I2,J,2)
			    Endif
              Enddo
		    Endif
	    Enddo
	  
	    Do I3 = 1,IJM
		    If(ISNBR(I3).EQ.2) Then
		    Do J = 1,CELL_POLYGEN(I3)
			    If((CELL_SIDE(I3,J,2).GT.0).AND.
     &               (ISNBR(CELL_SIDE(I3,J,2)).EQ.0)) Then
				    ISNBR(CELL_SIDE(I3,J,2)) = 3
				    NUM = NUM + 1
                      NBR(I,NUM) = CELL_SIDE(I3,J,2)
			    Endif
              Enddo
		    Endif
	    Enddo
	   	 
	    Do N = 1,IJM
		    ISNBR(N) = 0
          Enddo  
      Enddo
      
1000  Format (80A1)
      
      
      End Subroutine