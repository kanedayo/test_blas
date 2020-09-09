% cf. https://www.nag-j.co.jp/lapack/zposvx.htm
% zposvx : ~xは、エキスパートドライバ。様々な補助情報が得られえるが、メモリ使用量２倍になる。
A =[
 ( 3.23 +0.00i) ( 1.51 -1.92i) ( 1.90 +0.84i) ( 0.42 +2.50i)
 ( 1.51 +1.92i) ( 3.58 +0.00i) (-0.23 +1.11i) (-1.18 +1.37i)
 ( 1.90 -0.84i) (-0.23 -1.11i) ( 4.09 +0.00i) ( 2.33 -0.14i)
 ( 0.42 -2.50i) (-1.18 -1.37i) ( 2.33 +0.14i) ( 4.29 +0.00i)
 ];


B=[
 ( 3.93  -6.14i) ( 1.48  +6.58i)
 ( 6.17  +9.42i) ( 4.65  -4.75i)
 (-7.17 -21.83i) (-4.91  +2.29i)
 ( 1.99 -14.38i) ( 7.64 -10.79i)
 ];


X = A\B


%{
 ZPOSVX Example Program Results

 Solution(s)
                    1                 2
 1  ( 1.0000,-1.0000) (-1.0000, 2.0000)
 2  (-0.0000, 3.0000) ( 3.0000,-4.0000)
 3  (-4.0000,-5.0000) (-2.0000, 3.0000)
 4  ( 2.0000, 1.0000) ( 4.0000,-5.0000)

 Backward errors (machine-dependent)
       1.1E-16    5.3E-17

 Estimated forward error bounds (machine-dependent)
       6.0E-14    7.2E-14

 Estimate of reciprocal condition number
       6.6E-03

 A has not been equilibrated
 %}
