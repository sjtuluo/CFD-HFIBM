C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (28/07/2009)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c                                                                           #
c               Solving the trangle equations                               #
c                                                                           #
c############################################################################

	Subroutine ATRDE(B,N,M,D,L)
	Dimension B(M),D(N)
	Double Precision B,D
      
      
	L = 1
      
	If (M.NE.(3*N-2)) Then
          L = -1
	    Write(*,10)
	    Return
      Endif
      
10    Format(1X,'  ERR  ')
      
	Do 20 K = 1,N-1
          J = 3 * K - 2
	    If (ABS(B(J))+1.0.EQ.1.0) Then
	        L = 0
	        Write(*,10)
	        Return
	    Endif
	    B(J+1) = B(J+1)/B(J)
	    D(K) = D(K)/B(J)
	    B(J+3) = B(J+3)-B(J+2)*B(J+1)
	    D(K+1) = D(K+1)-B(J+2)*D(K)
20    Continue
      
	If (ABS(B(3*N-2))+1.0.EQ.1.0) Then
          L=0
	    Write(*,10)
	    Return
	Endif
	D(N) = D(N)/B(3*N-2)
      
	Do 30 K = N-1,1,-1
          D(K)=D(K)-B(3*K-1)*D(K+1)
30    Continue
      
	Return
      
	End
