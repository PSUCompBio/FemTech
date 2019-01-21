$# LS-DYNA Keyword file created by LS-PrePost(R) V4.5.23 - 17Aug2018
$# Created on Jan-17-2019 (13:57:46)
*KEYWORD  
*PART
$#                                                                         title
part1
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         1         1         1         0         0         0         0         0
*SECTION_SHELL_TITLE
sec1
$#   secid    elform      shrf       nip     propt   qr/irid     icomp     setyp
         1        12       1.0         2       1.0         0         0         1
$#      t1        t2        t3        t4      nloc     marea      idof    edgset
       0.0       0.0       0.0       0.0       0.0       0.0       0.0         0
*MAT_ELASTIC_TITLE
mat1
$#     mid        ro         e        pr        da        db  not used        
         1       1.0       1.0       1.0       0.0       0.0         0
*PART
$#                                                                         title
part2
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         2         1         2         0         0         0         0         0
*MAT_ELASTIC_TITLE
mat2
$#     mid        ro         e        pr        da        db  not used        
         2       1.0       1.0       0.2       0.0       0.0         0
*SET_SHELL_LIST_TITLE
set1
$#     sid       da1       da2       da3       da4   
         1       0.0       0.0       0.0       0.0
$#    eid1      eid2      eid3      eid4      eid5      eid6      eid7      eid8
         4         0         0         0         0         0         0         0
*SET_SHELL_LIST_TITLE
set2
$#     sid       da1       da2       da3       da4   
         2       0.0       0.0       0.0       0.0
$#    eid1      eid2      eid3      eid4      eid5      eid6      eid7      eid8
         2         1         3         0         0         0         0         0
*ELEMENT_SHELL
$#   eid     pid      n1      n2      n3      n4      n5      n6      n7      n8
       1       1       2       3       9       8       0       0       0       0
       2       2       6       8       9       7       0       0       0       0
       3       2       3       1       5       9       0       0       0       0
       4       2       5       4       7       9       0       0       0       0
*NODE
$#   nid               x               y               z      tc      rc  
       1             0.0             0.0             0.0       0       0
       2             0.0             1.0             0.0       0       0
       3             0.0             0.5             0.0       0       0
       4             1.0             0.0             0.0       0       0
       5             0.5             0.0             0.0       0       0
       6             1.0             1.0             0.0       0       0
       7             1.0             0.5             0.0       0       0
       8             0.5             1.0             0.0       0       0
       9             0.5             0.5             0.0       0       0
*END
