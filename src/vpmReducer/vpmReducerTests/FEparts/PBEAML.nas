BEGIN BULK
GRID       16352       0     1.5     3.8    2.66       0
GRID       21534       0     1.5     3.83.173333       0
GRID       21535       0     1.5     3.83.686667       0
GRID       16355       0     1.5     3.8     4.2       0
CBEAM      19721      18   16352   21534      0.     -1.      0.        +
+                       -2.727-3.0733612      0.-2.727-3.0733612      0.
CBEAM      19722      18   21534   21535      0.     -1.      0.        +
+                       -2.727-3.0733612      0.-2.727-3.0733612      0.
CBEAM      19723      18   21535   16355      0.     -1.      0.        +
+                       -2.727-3.0733612      0.-2.727-3.0733612      0.
PBEAML        18       1 MSCBML0       L                                +
+           .026     .08    .007    .006      0.
MAT1           1 2.06+117.923+10      .3   7850.      0.      0.
FORCE,1,16355,,10000.0,,1.0
SPC,,16352,123456,0.0
ENDDATA