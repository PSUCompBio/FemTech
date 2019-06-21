$# LS-DYNA Keyword file created by LS-PrePost(R) V4.5.23 - 17Aug2018
$# Created on Dec-09-2018 (10:45:38)
*KEYWORD  
*PART
$#                                                                         title
boxsolid
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         1         0         0         0         0         0         0         0
*ELEMENT_SOLID
$#   eid     pid      n1      n2      n3      n4      n5      n6      n7      n8
       1       1       1       2       4       3       5       6       8       7
       2       1       5       6       8       7       9      10      12      11
*NODE
$#   nid               x               y               z      tc      rc  
       1             0.0             0.0             0.0       0       0
       2             0.005           0.0             0.0       0       0
       3             0.0             0.005           0.0       0       0
       4             0.005           0.005           0.0       0       0
       5             0.0             0.0             0.0025    0       0
       6             0.005           0.0             0.0025    0       0
       7             0.0             0.005           0.0025    0       0
       8             0.005           0.005           0.0025    0       0
       9             0.0             0.0             0.005     0       0
      10             0.005           0.0             0.005     0       0
      11             0.0             0.005           0.005     0       0
      12             0.005           0.005           0.005     0       0
*END
