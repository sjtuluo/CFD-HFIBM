C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (21/02/2023)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                                                                           #     
C                    Save result data of IBM module                         #                                           
C                                                                           #
C############################################################################
      Subroutine IBMARCHIVE

      Include './Include/OCERM_INF'
      
      Double precision Cp1,Cp2,FD,FL,CD,CL
      Double precision UREF,ELFREF,PNREF
      Integer SAVE_ID(IJP,KB)
      Integer IJK,IJCELL
      Character*10 FN
      
      
C===========================================================================C      
C         Save results of time averaged
C===========================================================================C
      
      Do I = 1,IJP
          SAE(I) = SAE(I) + ELFV(I)
          Do K = 1,KB
              SAU(I,K) = SAU(I,K) + URV(I,K)
              SAV(I,K) = SAV(I,K) + VRV(I,K)
              SAW(I,K) = SAW(I,K) + WRV(I,K)
              SAP(I,K) = SAP(I,K) + PNV(I,K)
          Enddo
      Enddo

      If (NSTEP.EQ.IEND) Then
          Do I = 1,IJP
              ELFV(I) = SAE(I) / (IEND - ISTART)
              Do K = 1,KB
	            URV(I,K) = SAU(I,K) / (IEND - ISTART)
 	            VRV(I,K) = SAV(I,K) / (IEND - ISTART)
	            WRV(I,K) = SAW(I,K) / (IEND - ISTART)
                  PNV(I,K) = SAP(I,K) / (IEND - ISTART)
              Enddo
          Enddo
          
          Write(FN,'(F10.4)')  TIME * 24. + HOURSTAR
          IJK = IJP * KB
          IJCELL = IJM * KBM
          WRITE(ITF,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20,A20)')
     &        'ZONE T=', '"',FN,'"',' I=',IJK,' K=',IJCELL,',F=FEPOINT',
     &        ', ET=BRICK',', D=(FECONNECT)'
          Do K = 1, KB
          Do I = 1, IJP
              ZP = HP(I)+ELFV(I)
              If(ZP .LE. 0.0) Then
                  ZP = 0.001
                  ZZZ = Z(K) * ZP + (ZP - HP(I))
              Else
                  ZZZ = Z(K)*(HP(I)+ELFV(I))+ELFV(I)
              Endif	       
              Write(ITF,'(20E16.5)') 
     &             PXY(I,1), PXY(I,2), ZZZ, ELFV(I), PNV(I,K),
     &             URV(I,K), VRV(I,K), WRV(I,K), VORTICITYV(I,K),
     &             VORXV(I,K), VORYV(I,K), VORZV(I,K),OMIGAV(I,K), 
     &             QCRIV(I,K),KMV(I,K),TKEV(I,K),TDISSV(I,K)
          Enddo
          Enddo    
      Endif
      
C===========================================================================C      
C         Calculate and save Cp with THETA 
C===========================================================================C
      
      Do I = 1,IJM
          SAE2(I) = SAE2(I) + ELF(I)
          Do K = 1,KBM
              SAP2(I,K) = SAP2(I,K) + PN(I,K)
          Enddo
      Enddo
              
      If (NSTEP.EQ.IEND) Then
          Do I = 1,IJM
              SAE2(I) = SAE2(I) / (IEND - ISTART)
              Do K = 1,KBM
	            SAP2(I,K) = SAP2(I,K) / (IEND - ISTART)
              Enddo
          Enddo
          
          Do K = 1,KBM
          Do I = 1,IJM
          If (IDC(I,K).EQ.2) Then
              Cp1 = SAE2(I)
              Cp2 = SAP2(I,K)
              Write(IBMCP,'(5F16.8)') CXY(I,1), CXY(I,2), ZZ(k),
     &                                Cp1, Cp2
          Endif
          Enddo
          Enddo
      Endif

C===========================================================================C      
C         Calculate and save CD,CL with time 
C===========================================================================C

      If ((ISKILL.NE.0).AND.(Mod(NSTEP,ISKILL).EQ.0) 
     &	 .AND.(NSTEP.GT.ISTART)) Then
          TMIDDLE = TIME - DTI * DAYI / SKILLI / 2.
          TSAVE = (TMIDDLE * 24. + HOURSTAR) * 3600.
          FD = 0
          FL = 0
          Do I = 1,IJM
              If (IDC(I,K).GT.1.9) Then
                  FD = FD - BFX(I,K) * AREA(I)/(HC(I)+ELF(I))
                  FL = FL - BFY(I,K) * AREA(I)/(HC(I)+ELF(I))
              Endif
          Enddo
          CD = 20*FD/(UINF**2)
          CL = 20*FL/(UINF**2)
          
          Write(IBMCDCL,'(4F16.10)') TSAVE,CD,CL,UINF
      Endif     
          
C===========================================================================C      
C         Save results of ghost cells
C===========================================================================C
          
      Do K = 1,KB
      Do I = 1,IJP
          SAVE_ID(I,K) = 0
      Enddo
      Enddo
      
      Do K = 1,KBM
      Do I = 1,IJM
      If (IDC(I,K).GT.0.5) Then
          Do J1 = 1,CELL_POLYGEN(I)
          SAVE_ID(CELL_SIDE(I,J1,3),K) = 1
          SAVE_ID(CELL_SIDE(I,J1,3),K+1) = 1
          I2 = CELL_SIDE(I,J1,2)
          If (I2.GT.0) Then
              Do J2 = 1,CELL_POLYGEN(I2)
              SAVE_ID(CELL_SIDE(I2,J2,3),K) = 1
              SAVE_ID(CELL_SIDE(I2,J2,3),K+1) = 1
              Enddo
          Endif
          Enddo
      Endif
      Enddo
      Enddo
      
      Do K = 1,KB
          Do I = 1,IJP
              If (SAVE_ID(I,K).LT.0.5) Then
                  URV(I,K) = 0
                  VRV(I,K) = 0
                  WRV(I,K) = 0
                  ELFV(I) = 0
                  PNV(I,K) = 0  
                  TKEV(I,K) = 0
                  TDISSV(I,K) = 0
              Endif
          Enddo
      Enddo

      If((Mod(NSTEP,ISAVE).EQ.0).OR.(NSTEP.EQ.ISTART)) Then
          If(NSTEP.EQ.ISTART) Then
              Write(FN,'(F10.3)')  TIME * 24. + HOURSTAR
              IJK = IJP * KB
              IJCELL = IJM * KBM
              Write(IBMDATA,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20)')
     &        'ZONE T=', '"',FN,'"',' I=',IJK,' K=',IJCELL,',F=FEPOINT',
     &        ', ET=BRICK'
     
              Do K = 1, KB
              Do I = 1, IJP
                  ZP = HP(I)+ELFV(I)
                  If(ZP .LE. 0.0) Then
                      ZP = 0.001
                      ZZZ = Z(K) * ZP + (ZP - HP(I))
                  Else
                      ZZZ = Z(K)*(HP(I)+ELFV(I))+ELFV(I)
                  Endif
                  Write(IBMDATA,'(10E16.8)') 
     &                PXY(I,1), PXY(I,2), ZZZ, ELFV(I), PNV(I,K),
     &                URV(I,K), VRV(I,K), WRV(I,K),TKEV(I,K),TDISSV(I,K)
              Enddo
              Enddo
              
              Do K = 1, KBM 
              Do I = 1, IJM
                  Write(IBMDATA,'(8I8)') 
     &                IJP * (K-1) + CELL_INDEX(I,1), 
     &                IJP * (K-1) + CELL_INDEX(I,2), 
     &                IJP * (K-1) + CELL_INDEX(I,3), 
     &                IJP * (K-1) + CELL_INDEX(I,4), 
     &                IJP * (K) + CELL_INDEX(I,1), 
     &                IJP * (K) + CELL_INDEX(I,2), 
     &                IJP * (K) + CELL_INDEX(I,3), 
     &                IJP * (K) + CELL_INDEX(I,4)
              Enddo   
              Enddo
	    Else if (NSTEP.GE.5000) Then
	        Write(FN,'(F10.4)')  TIME * 24. + HOURSTAR
		    IJK = IJP * KB
	        IJCELL = IJM * KBM
	        Write(IBMDATA,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20,A20)')
     &        'ZONE T=', '"',FN,'"',' I=',IJK,' K=',IJCELL,',F=FEPOINT',
     &        ', ET=BRICK',', D=(FECONNECT)'
     
	        Do K = 1, KB
              Do I = 1, IJP
                  ZP = HP(I)+ELFV(I)
                  If(ZP .LE. 0.0) Then
                      ZP = 0.001
                      ZZZ = Z(K) * ZP + (ZP - HP(I))
                  Else
                      ZZZ = Z(K)*(HP(I)+ELFV(I))+ELFV(I)
                  Endif	       
                  Write(IBMDATA,'(10E16.8)') 
     &                PXY(I,1), PXY(I,2), ZZZ, ELFV(I), PNV(I,K),
     &                URV(I,K), VRV(I,K), WRV(I,K),TKEV(I,K),TDISSV(I,K)
              Enddo
	        Enddo
	    Endif
      Endif
      
      End Subroutine