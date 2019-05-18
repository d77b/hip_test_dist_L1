Fortran Code for L-p Distance Statistic
The following two subroutines may be used to calculate the L-1 and L-2 distances between two density estimates. The L-1 statistic has been used very successfully in testing the hypothesis of equality of distributions from two samples. The reader is refered to the manuscript, "Hypothesis Testing Using an L-1 Distance Bootstrap" by David Allen, Eastern Oregon State College, appearing in the May, 1997 issue of The American Statistician.

The source code (FORTRAN language) for each of the two subroutines may be copied by clicking on the file name below, then choosing File...Save As... from your browser menu. The code may also be viewed below.

LPDIST.FOR
Code for calculating the L1 and L2 distances between two density estimates
BUILDENS.FOR
Code for building a density estimate from an ordered list of a sample
General Introduction to the Routines
For a given sample, the data must first be ordered. There are numerous sorting procedures available and one is not provided below.	From the ordered sample, the subroutine BUILDENS will build the density estimate as described in the article "Hypothesis Testing Using an L-1 Distance Bootstrap" by David Allen, Eastern Oregon State College, appearing in the May, 1997 issue of The American Statistician. The density estimate is a step function, specified by a sequence of pairs (S(i),F(i)) stored in the arrays S and F.	These pairs identify the left endpoints of each "step", that is, f(x)=F(i) for x in [S(i),S(i+1)). The subroutine LPDIST calculates the L-1, L-2, and L-1/2 distances between two density estimates created by BUILDENS. The code may be copied from this web site and may easily be modified to adapt to other problems of interest. While the code has been carefully tested, we cannot guarantee against unforeseen problems.

Source Code



      SUBROUTINE BUILDENS ( X, N, C, S, F, NF)

C  This subroutrine builds a Kernel density estimate from an ordered
C  set of observations X(1) ... X(N).  We use the Kernel function :
C  K(x) = 1/2  ( = 0 ) for ABS(x) < 1  (ABS(x) > 1).
C  With this kernel, the density estimate is a step function of the form
C
C			       ___ N
C	      f(x) =  (C/n) *  \ 	K (C*(x-X(i))	.
C			       /__ i=1
C
C  Inputs to the subroutine are the observation array X, its length N, and
C  the parameter C associated with the kernel estimate.  A recommended
C  value for C is (N**0.2)/(2.0*SD), where SD is the sample standard deviation.
C
C  Outputs from the subroutine are the arrays S and F of length NF which
C  define the density estimate.  Note that each observation generates two
C  discontinuities in the density estimate (one step up and one step down)
C  so the array length NF will generally satisfy  NF = 2 + 2*N.  THIS FACT
C  IS IMPORTANT IN ARRAY DECLARATION.	 The estimate
C  is defined in these arrays by  f(x) = F(i) for x in [ S(i) , S(i+1) ),
C  i ranging fro 1 to NF.  As a convention, F(1) = F(NF) = F(NF-1) = 0.
C  This permits some manipulation of the quantities S(1) and S(NF) which
C  allows adjustment of the "support" of the density estimate.	In
C  particular, S(1) may be assigned any value less than S(2) and S(NF) may
C  be assigned any value greater than S(NF-1) without changing the density.
C  Such adjustments permit the supports of two density estimates to coincide.
C
      IMPLICIT NONE
      REAL*4 X(1),    ! input array of ordered observations
     & C,	      ! input parameter defining the kernel width
     & HSIZE,	      ! vertical step size in the estimate (equals 'C/2N')
     & R,	      ! radius of the Kernel function  (equals	'1/C')
     & S(1),	      ! position along the x-axis
     & F(1),	      ! heights of the function
     & CP	      ! current position

      INTEGER N,      ! the length of X(i), the array of observations
     & NF,	      ! the length of arrays S(i) and F(i) (generally NF=2+2*N)
     & ICP,	      ! index for the current position
     & IL,	      ! index for the lower limit of the Kernel domain
     & IH	      ! index for the upper limit of the Kernel domain

C  In the initialization, it is necessary to create a terminal value X(N+1)
C  in the array X(i).  This is used in the check for completion(?!!).  The
C  algorithm steps through the locations for jump discontinuities in the
C  density estimate and determines the value the density estimate takes on
C  at that jump.  Note the jumps occur only at X(i)+R or X(i)-R.
C  The support of the density estimate begins at S(2) = X(1) - R but we
C  define the density estimate on a wider support, beginning at S(1).
C  It is convenient later for the first interval to have function value 0.

      HSIZE = C / FLOAT(2*N)
      R = 1.0 / C
      X(N+1) = X(N) + 4*R	    ! note X(N+1)-R < X(N)+R
      S(2) = X(1) - R
      S(1) = S(2) - 10. 	    !  the choice "-10" is arbitrary
      F(1) = 0.
      F(2) = HSIZE
      ICP = 2
      CP = S(2)
      IL = 2
      IH = 1

C  Here we begin the loop.  Note that at the last pass, the last discontinuity
C  has been established at X(N)+R.  Thus the high or upper endpoint has been
C  processed for the index N so the index IH will then be set to N+1.  For this
C  reason, the check for completion is that IH exceeds N.

      DO WHILE (IH.LE.N)
	 IF ( X(IL)-R .LE. X(IH)+R ) THEN
	    IF ( CP .EQ. X(IL)-R ) THEN 	!  there is no new jump point
	       F(ICP) = F(ICP) + HSIZE		!  but the current function
	    ELSE				!  value must be increased
	       ICP = ICP + 1
	       CP = X(IL) - R			!  set the new jump point and
	       S(ICP) = CP			!  the value of the function
	       F(ICP) = F(ICP-1) + HSIZE	!  at that jump
	    ENDIF
	    IL = IL + 1
	 ELSE
	    IF ( CP .EQ. X(IH)+R ) THEN
	       F(ICP) = F(ICP) - HSIZE		!  this code is analogous to
	    ELSE				!  the preceding code except
	       ICP = ICP + 1			!  the step function is decreasing
	       CP = X(IH) + R			!  rather than increasing at these
	       S(ICP) = CP			!  points of discontinuity
	       F(ICP) = F(ICP-1) - HSIZE
	    ENDIF
	    IH = IH + 1
	 ENDIF
      END DO
      NF = ICP + 1				!  for the purpose of aligning
      S(NF) = S(ICP) + 10.			!  the interval of support, the
      F(NF) = 0.				!  last element of S(i) is set
      RETURN					!  to any value satisfying
      END					!  S(ICP) < S(NF)








      SUBROUTINE LPDIST ( S, F, NF, T, G, NG, L1, L2, LH )
C
C  This subroutine computes the Lp distances between two functions F and G.
C  The choices for p are 1, 2, and 1/2 yielding measures L1, L2, and LH.  The
C  function F is given by the two arrays S(i) and F(i) of length NF as
C  produced by the kernel density estimation subroutine  BUILDENS.  Similarly
C  G is given by T(j) and G(j) of length NG.  Each function is a step function
C  so the L! distance is a sum of rectangular areas.  The algorithm steps
C  along the x-axis identifying the lower and upper limits (LOW and HIGH) of
C  each such rectangle and the corresponding contribution to the L1 distance.
C
      IMPLICIT NONE
      INTEGER NF,	 !  length of the arrays S and F
     & NG,		 !  length of the arrays T and G
     & INEXTF,		 !  index used to point to the next step of F
     & INEXTG		 !  index used to point to the next step of G

      REAL*4 S(1),	 !  array of x-values for step locations in F
     & F(1),		 !  array of function heights after steps in F
     & T(1),		 !  array of x-values for step locations in G
     & G(1),		 !  array of function heights after steps in G
     & FCUR,		 !  variable holding the current value of F
     & GCUR,		 !  variable holding the current value of G
     & LOW,		 !  the lower limit of the interval being integrated
     & HIGH,		 !  the upper limit ...
     & END,		 !  the far right endpoint of integration
     & L1,		 !  the value of the distance function
     & L2,		 !  the value of the distance function
     & LH		 !  the value of the distance function


      S(1) = MIN( S(1) , T(1) )
      T(1) = S(1)
      S(NF) = MAX( S(NF-1) , T(NG-1) )
      T(NG) = S(NF)
      END = S(NF)
      FCUR = F(1)
      GCUR = G(1)
      INEXTF = 2
      INEXTG = 2
      LOW = S(1)
      L1 = 0.0
      L2 = 0.0
      LH = 0.0
      DO WHILE ( LOW .LT. END )
	 IF (T(INEXTG) .LT. S(INEXTF)) THEN
	    HIGH = T(INEXTG)
	    L1 = L1 + ABS(FCUR-GCUR)*(HIGH-LOW)
	    L2 = L2 + ABS(FCUR-GCUR)*ABS(FCUR-GCUR)*(HIGH-LOW)
	    LH = LH + SQRT(ABS(FCUR-GCUR))*(HIGH-LOW)
	    GCUR = G(INEXTG)
	    INEXTG = INEXTG + 1
	 ELSE
	    HIGH = S(INEXTF)
	    L1 = L1 + ABS(FCUR-GCUR)*(HIGH-LOW)
	    L2 = L2 + ABS(FCUR-GCUR)*ABS(FCUR-GCUR)*(HIGH-LOW)
	    LH = LH + SQRT(ABS(FCUR-GCUR))*(HIGH-LOW)
	    FCUR = F(INEXTF)
	    INEXTF = INEXTF + 1
	 ENDIF
	 LOW = HIGH
      END DO
      L2=SQRT(L2)
      LH=LH*LH
      RETURN
      END




Go to Allen's Home Page