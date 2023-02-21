C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (21/02/2023)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                                                                           #     
C                 Main module of Immersed Boundary Method                   #                                           
C                                                                           #
C############################################################################
      Subroutine IBM
      
      Include './Include/OCERM_INF'
      
      
      If (NSTEP.EQ.ISTART) Then
          Write(*,1000) IBMDOMAIN,
     &    'Ray-casting method (2-2D Boundary,3-3D Boundary)'
          Write(*,1000) IBMSCHEME,
     &    'Boundary mode (0-Static Boundary,1-Moving Boundary)'
          Write(*,1000) IBGCMODE,
     &    'Ghost cell mode (0-Without pressure BC,1-Puessrue BC)'
          Write(*,1000) IBFSCHEME,
     &    'Force term mode (1-Ghost Cell MODE,2-Direct Force MODE)'
      Endif

C===========================================================================C      
C         Different index of cells detecting
C===========================================================================C

      If (IBMSCHEME.GE.0) Then
          Call IBMALIDC                         
          Call IBMALIDP
      Endif
      
C===========================================================================C      
C         Ghost cell & Direct force immersed boundary method  
C===========================================================================C
      
      If (IBFSCHEME.EQ.1) Call IBMALGC
      If (IBFSCHEME.EQ.2) Call IBMALDF
      
C===========================================================================C      
C         Visualization of immersed boundary method 
C===========================================================================C

      Call IBMVISUAL(1,100) 
      
      If (IBMSCHEME.EQ.0) Then
          Write(*,*) 'Static boundary mode.'
          IBMSCHEME = -1
      Else if (IBMSCHEME.EQ.1) Then
          Write(*,*) 'Moving boundary mode.'
          IBMSCHEME = 2
      Endif
         
C===========================================================================C     
C         Update moving boundary 
C===========================================================================C
     
      If (IBMSCHEME.GE.1) Call IBMRGD
      
      If (NSTEP.EQ.ISTART) Pause
      
      
1000  Format (I3,2X,A)
      
      
      End Subroutine