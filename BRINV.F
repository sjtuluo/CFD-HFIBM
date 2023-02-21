C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (28/07/2009)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c                                                                           #
c               inverse of the coefficient matrix                           #
c                                                                           #
c############################################################################
      
	Subroutine BRINV(A,N,L)
	Dimension A(N,N),IS(N),JS(N)
	Double Precision A,T,D
      
      
	L = 1
      
	Do 100 K = 1,N
          D = 0.0
          
	    Do 10 I = K,N
	    Do 10 J = K,N
              If (Abs(A(I,J)).GT.D) Then
                  D = Abs(A(I,J))
	            IS(K) = I
	            JS(K) = J
	        Endif
10        Continue
          
	    If (D+1.0.EQ.1.0) Then
	        L = 0
	        Write(*,20)
	        Return
          Endif
          
20        Format(1X,'ERR**NOT INV')
          
	    Do 30 J = 1,N
	        T = A(K,J)
	        A(K,J) = A(IS(K),J)
	        A(IS(K),J) = T
30        Continue
          
	    Do 40 I = 1,N
	        T = A(I,K)
	        A(I,K) = A(I,JS(K))
	        A(I,JS(K)) = T
40        Continue
          
	    A(K,K) = 1/A(K,K)
	    Do 50 J = 1,N
	        If (J.NE.K) Then
                  A(K,J) = A(K,J)*A(K,K)
	        Endif
50        Continue
          
	    Do 70 I = 1,N
	        If (I.NE.K) Then
	        Do 60 J = 1,N
	            If (J.NE.K) Then
	                A(I,J) = A(I,J)-A(I,K)*A(K,J)
	            Endif
60	        Continue
	        Endif
70        Continue
          
	    Do 80 I = 1,N
	        If (I.NE.K) Then
	            A(I,K) = -A(I,K)*A(K,K)
	        Endif
80        Continue
          
100   Continue
      
	Do 130 K = N,1,-1
          Do 110 J = 1,N
	        T = A(K,J)
	        A(K,J) = A(JS(K),J)
	        A(JS(K),J) = T
110       Continue
          
	    Do 120 I = 1,N
	        T = A(I,K)
	        A(I,K) = A(I,IS(K))
	        A(I,IS(K)) = T
120       Continue
          
130   Continue
      
	Return
      
	End
