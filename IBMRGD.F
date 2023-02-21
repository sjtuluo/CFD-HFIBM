C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (21/02/2023)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c                                                                           #
c        Rigid body dynamics & update coordinates of moving boundary        #                                           
c                                                                           #
c############################################################################
      Subroutine IBMRGD

      Include './Include/OCERM_INF'
      
      Type QT
          Real s, x, y, z
      End Type

      Real(Kind=8) MAT_ROT(3,3),MAT_ROT_T(3,3),MAT_I(3,3),MAT_I1(3,3)
      Real(Kind=8) MAT_IR(3,3),VEC_TOR(3),RI(3),RRI(3),RRN(3),PLANE(4)
      Real(Kind=8) MAT_TEMP(3,3),MAT_TEMP2(3,3),VEC_TEMP(3),VEC_TEMP2(3)
      Real(Kind=8) VEC_J(3),X_SUM,Y_SUM,Z_SUM,QT_MAG,MI,UDE,WDE
      Integer NDE

      
      Type(QT) :: QT0, QT1, QTW

      Do I = 1,3
          VEC_TOR(I) = 0
      Enddo
      
      
      NDE = 10
      If (Mod(NSTEP,NDE) .EQ. 0.0) Then
          UDE = 0.999
          WDE = 0.999
      Else
          UDE = 1
          WDE = 1
      Endif
      
      MI = 0.1
      If (NSTEP.EQ.ISTART) Then
      Do NB = 1,NUMIB
      FRN(NB) = 0.9
      RGF(NB,1) = 0
      RGF(NB,2) = 0
      RGF(NB,3) = -9.8
      Do N = 1,IBN(NB)
          RI(1) = RI0(NB,N,1)
          RI(2) = RI0(NB,N,2)
          RI(3) = RI0(NB,N,3)
          RGIR(NB,1,1) = RGIR(NB,1,1) + MI*(RI(2)*RI(2)+RI(3)*RI(3))
          RGIR(NB,1,2) = RGIR(NB,1,2) - MI*RI(1)*RI(2)
          RGIR(NB,1,3) = RGIR(NB,1,3) - MI*RI(1)*RI(3)
          RGIR(NB,2,1) = RGIR(NB,2,1) - MI*RI(2)*RI(1)
          RGIR(NB,2,2) = RGIR(NB,2,2) + MI*(RI(1)*RI(1)+RI(3)*RI(3))
          RGIR(NB,2,3) = RGIR(NB,2,3) - MI*RI(2)*RI(3)
          RGIR(NB,3,1) = RGIR(NB,3,1) - MI*RI(3)*RI(1)
          RGIR(NB,3,2) = RGIR(NB,3,2) - MI*RI(3)*RI(2)
          RGIR(NB,3,3) = RGIR(NB,3,3) + MI*(RI(1)*RI(1)+RI(2)*RI(2)) 
      Enddo
      Enddo
      Endif
          
      Do 100 NB = 1,NUMIB 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C  
      Do I = 1,3
      Do J = 1,3
          MAT_IR(I,J) = RGIR(NB,I,J)
      Enddo
      Enddo
      
      RGM(NB) = MI * IBN(NB)
      
      QT_MAG = Sqrt(RGQT(NB,1)**2 + RGQT(NB,2)**2 + 
     &              RGQT(NB,3)**2 + RGQT(NB,4)**2)
     &              + 0.00000001            
      QT0.s = RGQT(NB,1) / QT_MAG 
      QT0.x = RGQT(NB,2) / QT_MAG 
      QT0.y = RGQT(NB,3) / QT_MAG 
      QT0.z = RGQT(NB,4) / QT_MAG 
      
C===========================================================================C     
C         Step 1: Update velocity and angle velocity
C===========================================================================C  

      ! Update velocity
      Do I = 1,3
          RGU(NB,I) = RGU(NB,I) + DTI * RGF(NB,I) / RGM(NB)
      Enddo
      
      MAT_ROT(1,1) = QT0.s * QT0.s + QT0.x * QT0.x -
     &               QT0.y * QT0.y - QT0.z * QT0.z
      MAT_ROT(1,2) = 2 * (QT0.x*QT0.y - QT0.s*QT0.z)
      MAT_ROT(1,3) = 2 * (QT0.x*QT0.z + QT0.s*QT0.y)
      MAT_ROT(2,1) = 2 * (QT0.x*QT0.y + QT0.s*QT0.z)
      MAT_ROT(2,2) = QT0.s * QT0.s - QT0.x * QT0.x +
     &               QT0.y * QT0.y - QT0.z * QT0.z
      MAT_ROT(2,3) = 2 * (QT0.y*QT0.z - QT0.s*QT0.x)
      MAT_ROT(3,1) = 2 * (QT0.x*QT0.z - QT0.s*QT0.y)
      MAT_ROT(3,2) = 2 * (QT0.y*QT0.z + QT0.s*QT0.x)
      MAT_ROT(3,3) = QT0.s * QT0.s - QT0.x * QT0.x -
     &               QT0.y * QT0.y + QT0.z * QT0.z
      
      MAT_ROT_T = TRANSPOSE(MAT_ROT) 
      
      MAT_TEMP = MATMUL(MAT_ROT,MAT_IR)
      MAT_I = MATMUL(MAT_TEMP,MAT_ROT_T)
      
      ! Update angle velocity
      Call MAT_INV(MAT_I,MAT_I1)
      Call MUL_VEC(MAT_I1,VEC_TOR,VEC_TEMP) 
      
      Do I = 1,3
          RGW(NB,I) = RGW(NB,I) + DTI * VEC_TEMP(I)
      Enddo

C===========================================================================C     
C        Step 2: Collision Detection
C===========================================================================C  

      PLANE(1) = 0
      PLANE(2) = 0
      PLANE(3) = 1
      PLANE(4) = -0.6

      Call IBMCLBP(NB,MAT_ROT,MAT_I1,PLANE)
      
C===========================================================================C     
C         Step 3: Update position and quaternion
C===========================================================================C 
      ! Update position of gravity center      
      Do I = 1,3
          RGX(NB,I) = RGX(NB,I) + DTI * RGU(NB,I)
          RGU(NB,I) = UDE * RGU(NB,I)
          RGW(NB,I) = WDE * RGW(NB,I)
      Enddo      
      
      ! Update w and q
      QTW.s = 0
      QTW.x = 0.5 * DTI * RGW(NB,1)
      QTW.y = 0.5 * DTI * RGW(NB,2)
      QTW.z = 0.5 * DTI * RGW(NB,3)

      Call MUL_QT(QTW,QT0,QT1)
      QT1.s = QT1.s + QT0.s
      QT1.x = QT1.x + QT0.x
      QT1.y = QT1.y + QT0.y
      QT1.z = QT1.z + QT0.z
      QT_MAG = Sqrt(QT1.s**2 + QT1.x**2 + QT1.y**2 + QT1.z**2)
     &              + 0.00000001
      RGQT(NB,1) = QT1.s / QT_MAG
      RGQT(NB,2) = QT1.x / QT_MAG
      RGQT(NB,3) = QT1.y / QT_MAG
      RGQT(NB,4) = QT1.z / QT_MAG
      
      ! Update position of vertices and normal vector
      Do N = 1,IBN(NB)
          RI(1) = RI0(NB,N,1)
          RI(2) = RI0(NB,N,2)
          RI(3) = RI0(NB,N,3)
          Call MUL_VEC(MAT_ROT,RI,RRI)
          BXY(NB,N,1) = RGX(NB,1) + RRI(1)
          BXY(NB,N,2) = RGX(NB,2) + RRI(2)
          BXY(NB,N,3) = RGX(NB,3) + RRI(3)
          RI(1) = RN0(NB,N,1)
          RI(2) = RN0(NB,N,2)
          RI(3) = RN0(NB,N,3)
          Call MUL_VEC(MAT_ROT,RI,RRN)
          BXN(NB,N,1) = RRN(1)
          BXN(NB,N,2) = RRN(2)
          BXN(NB,N,3) = RRN(3)
      Enddo
      
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
100   Continue

      End Subroutine
      
      
C===========================================================================C
C      
C         Operation of Matrix and Vector
C      
C===========================================================================C  
      
C---------------------------------------------------------------------------C     
C         Matrix multiply vector
C---------------------------------------------------------------------------C  
      
      Subroutine MUL_VEC(MAT_A,VEC_B,VEC_C)
      
      Real(kind=8) MAT_A(3,3), VEC_B(3), VEC_C(3)
      Do I = 1,3
          VEC_C(I) = MAT_A(I,1)*VEC_B(1) + MAT_A(I,2)*VEC_B(2)
     &               + MAT_A(I,3)*VEC_B(3)
      Enddo
      
      End Subroutine
      
C---------------------------------------------------------------------------C     
C         Dot product of vectors
C---------------------------------------------------------------------------C  
      
      Subroutine DOT_VEC(VEC_A,VEC_B,SCA_C)
      
      Real(kind=8) VEC_A(3), VEC_B(3), SCA_C
      SCA_C = VEC_A(1) * VEC_B(1) + VEC_A(2) * VEC_B(2) +
     &        VEC_A(3) * VEC_B(3)
      
      End Subroutine   
      
C---------------------------------------------------------------------------C     
C         Cross product of vectors
C---------------------------------------------------------------------------C  
      
      Subroutine CRP_VEC(VEC_A,VEC_B,VEC_C)
      
      Real(kind=8) VEC_A(3), VEC_B(3), VEC_C(3)
      VEC_C(1) = VEC_A(2)*VEC_B(3) - VEC_B(2)*VEC_A(3)
      VEC_C(2) = VEC_A(3)*VEC_B(1) - VEC_B(3)*VEC_A(1)
      VEC_C(3) = VEC_A(1)*VEC_B(2) - VEC_B(1)*VEC_A(2)
      
      End Subroutine
      
C---------------------------------------------------------------------------C     
C         Multiplication of quaternions
C---------------------------------------------------------------------------C  

      Subroutine MUL_QT(QT_A,QT_B,QT_C)
      
      Type QT
          Real s, x, y, z
      End Type
      Type(QT) :: QT_A, QT_B, QT_C
      QT_C.s = QT_A.s * QT_B.s - QT_A.x * QT_B.x
     &          - QT_A.y * QT_B.y - QT_A.z * QT_B.z
      QT_C.x = QT_A.s * QT_B.x + QT_B.s * QT_A.x
     &          + QT_A.y * QT_B.z - QT_B.y * QT_A.z
      QT_C.y = QT_A.s * QT_B.y + QT_B.s * QT_A.y
     &          + QT_A.z * QT_B.x - QT_B.z * QT_A.x
      QT_C.z = QT_A.s * QT_B.z + QT_B.s * QT_A.z
     &          + QT_A.x * QT_B.y - QT_B.x * QT_A.y
      
      End Subroutine

C---------------------------------------------------------------------------C     
C         Multiplication of quaternions
C---------------------------------------------------------------------------C 

      Subroutine MAT_INV(MAT_A, MAT_B)
      
      Real(kind=8) MAT_A(3,3), MAT_B(3,3), MAT_A1(3,3)
      Real(kind=8) DET_A
      
      MAT_A1(1,1) = MAT_A(2,2) * MAT_A(3,3) - MAT_A(2,3) * MAT_A(3,2)
      MAT_A1(1,2) = MAT_A(3,1) * MAT_A(2,3) - MAT_A(2,1) * MAT_A(3,3)
      MAT_A1(1,3) = MAT_A(2,1) * MAT_A(3,2) - MAT_A(3,1) * MAT_A(2,2)
      MAT_A1(2,1) = MAT_A(1,3) * MAT_A(3,2) - MAT_A(1,2) * MAT_A(3,3)
      MAT_A1(2,2) = MAT_A(1,1) * MAT_A(3,3) - MAT_A(1,3) * MAT_A(3,1)
      MAT_A1(2,3) = MAT_A(1,2) * MAT_A(3,1) - MAT_A(1,1) * MAT_A(3,2)
      MAT_A1(3,1) = MAT_A(1,2) * MAT_A(2,3) - MAT_A(1,3) * MAT_A(2,2)
      MAT_A1(3,2) = MAT_A(1,3) * MAT_A(2,2) - MAT_A(1,1) * MAT_A(2,3)
      MAT_A1(3,3) = MAT_A(1,1) * MAT_A(2,2) - MAT_A(1,2) * MAT_A(2,1)
      DET_A = MAT_A(1,1) * MAT_A(2,2) * MAT_A(3,3)
     &        + MAT_A(1,2) * MAT_A(2,3) * MAT_A(3,1)
     &        + MAT_A(1,3) * MAT_A(3,2) * MAT_A(2,1)
     &        - MAT_A(1,3) * MAT_A(2,2) * MAT_A(3,1)
     &        - MAT_A(1,2) * MAT_A(2,1) * MAT_A(3,3)
     &        - MAT_A(1,1) * MAT_A(3,2) * MAT_A(2,3)
      
      If (DET_A.EQ.0) DET_A = 0.000001
      Do I = 1,3
          Do J = 1,3
              MAT_B(I,J) = MAT_A1(I,J) / DET_A 
          Enddo
      Enddo
      
      End Subroutine
