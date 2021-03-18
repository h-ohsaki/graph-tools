use Test::Simple 'tests' : 3

use Graph::Enhanced

g = Graph(directed=True, multiedged=True)
buf = split(/\n/, <<EOF)
LEDA.GRAPH
void
void
16
|{}|
|{}|
|{}|
|{}|
|{}|
|{}|
|{}|
|{}|
|{}|
|{}|
|{}|
|{}|
|{}|
|{}|
|{}|
|{}|
30
1 5 0 |{}|
1 1 0 |{}|
1 12 0 |{}|
2 5 0 |{}|
2 4 0 |{}|
2 1 0 |{}|
2 10 0 |{}|
3 14 0 |{}|
3 9 0 |{}|
4 14 0 |{}|
5 6 0 |{}|
5 7 0 |{}|
6 15 0 |{}|
6 8 0 |{}|
6 8 0 |{}|
8 16 0 |{}|
8 16 0 |{}|
8 13 0 |{}|
9 6 0 |{}|
9 3 0 |{}|
9 15 0 |{}|
10 13 0 |{}|
11 16 0 |{}|
12 6 0 |{}|
13 9 0 |{}|
14 16 0 |{}|
14 2 0 |{}|
14 16 0 |{}|
15 10 0 |{}|
16 3 0 |{}|
# version string
Graph::EnhancedWin 1.400000
# scaling  wxmin  wymin  wxmax  wymax
0.8823242 -10 -1.867188 499.9854 513.7847
# node label font and size
0 13.42105
# edge label font and size
0 12.46241
# node index format
d
# edge index format
d
# multi-edge distance
4.793233
#
# node infos
# x y shape bclr(r,g,b) bwidth r1 r2 clr(r,g,b) ltype lclr(r,g,b) lpos lstr
405.1759 232.5947 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
443.3203 337.2901 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
311.6545 100.4818 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
140.5621 363.9527 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
366.314 54.79785 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
243.8891 244.7433 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
240.4881 369.7675 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
175.8479 247.6566 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
195.6626 374.3596 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
338.0766 457.1196 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
397.3704 306.3133 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
53.89083 100.5278 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
99.79392 365.232 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
46.66504 331.2502 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
165.4018 66.97121 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
84.05275 113.269 0 0 0 0 0.575188 13.42105 13.42105 255 255 228 4 -1 -1 -1 4
#
# edge infos
# width clr(r,g,b) shape style dir ltype lclr(r,g,b) lpos sanch tanch poly lstr
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (402.31,219.4832) (369.1799,67.90936)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 5 (405.1759,219.1736) (405.1759,205.7526) (378.3338,205.7526) (378.3338,232.5947) (391.7548,232.5947)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (392.6133,227.8717) (66.45341,105.2508)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (439.7906,324.3415) (369.8438,67.74643)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (429.951,338.4675) (153.9314,362.7753)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (438.726,324.6799) (409.7702,245.2049)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (434.4638,347.3741) (346.9332,447.0357)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (301.5333,109.2958) (56.78615,322.4362)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (304.2977,111.7069) (198.6056,361.2652)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (127.8878,359.5385) (59.3394,335.6645)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (359.0432,66.07878) (251.16,233.4623)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (361.3351,67.26119) (245.467,357.3042)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (238.4684,232.4656) (170.8224,79.24888)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (230.7984,247.7026) (189.1437,249.4862)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (230.5933,242.9138) (188.9386,244.6973)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (166.4206,238.1041) (89.52206,125.5251)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (170.3786,235.4006) (93.48006,122.8215)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (168.5585,258.9256) (107.0833,353.963)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (200.3427,361.781) (239.209,257.3219)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (203.0193,363.1345) (308.7115,113.5762)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (194.3477,361.0031) (166.7167,80.3277)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (325.5544,452.2907) (112.3162,370.0609)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (385.9441,299.2731) (95.4791,120.3091)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (64.58113,108.6421) (233.1988,236.629)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (113.1546,366.5041) (182.302,373.0875)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (46.53527,317.8298) (79.45827,125.8792)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (60.08454,331.4546) (429.9008,337.0857)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (51.25952,318.6401) (84.18251,126.6895)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (170.8336,79.24396) (332.6449,444.8469)
0.9586466 0 0 0 0 0 1 1 0 0 0 5 (0,0) (0,0) 2 (97.45267,112.5162) (298.2545,101.2346)
EOF
g.import_graph('gw', buf)
ok(len(g.vertices()) == 16)
ok(len(g.edges()) == 30)
ok(g.has_edge(8, 16))
