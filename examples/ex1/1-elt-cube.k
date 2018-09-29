$# LS-DYNA Keyword file created by LS-PrePost(R) V4.5.23 - 17Aug2018
$# Created on Sep-29-2018 (04:08:04)
*KEYWORD
*PART
$#                                                                         title
boxsolid
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         1         0         0         0         0         0         0         0
*ELEMENT_SOLID
$#   eid     pid      n1      n2      n3      n4      n5      n6      n7      n8
       1       1       1       2       5       4       7       8      11      10
       2       1       2       3       6       5       8       9      12      11
*NODE
$#   nid               x               y               z      tc      rc
       1             0.0             0.0             0.0       0       0
       2            50.0             0.0             0.0       0       0
       3           100.0             0.0             0.0       0       0
       4             0.0           100.0             0.0       0       0
       5            50.0           100.0             0.0       0       0
       6           100.0           100.0             0.0       0       0
       7             0.0             0.0           100.0       0       0
       8            50.0             0.0           100.0       0       0
       9           100.0             0.0           100.0       0       0
      10             0.0           100.0           100.0       0       0
      11            50.0           100.0           100.0       0       0
      12           100.0           100.0           100.0       0       0
*END
