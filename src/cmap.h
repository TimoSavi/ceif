/* A table for mapping ascii char values to unified values. 
 * lower case chars (a-x) are mapped to same value as upper case (A-Z).
 */

struct char_mapping
{
    unsigned char map_value;
};

struct char_mapping cmap[] =
{
{1},    // orig dec=0, char=NUL
{1},    // orig dec=1, char=SOH
{2},    // orig dec=2, char=STX
{3},    // orig dec=3, char=ETX
{4},    // orig dec=4, char=EOT
{5},    // orig dec=5, char=ENQ
{6},    // orig dec=6, char=ACK
{7},    // orig dec=7, char=BEL
{8},    // orig dec=8, char=BS
{9},    // orig dec=9, char=HT
{10},    // orig dec=10, char=LF
{11},    // orig dec=11, char=VT
{12},    // orig dec=12, char=FF
{13},    // orig dec=13, char=CR
{14},    // orig dec=14, char=SO
{15},    // orig dec=15, char=SI
{16},    // orig dec=16, char=DLE
{17},    // orig dec=17, char=DC1
{18},    // orig dec=18, char=DC2
{19},    // orig dec=19, char=DC3
{20},    // orig dec=20, char=DC4
{21},    // orig dec=21, char=NAK
{22},    // orig dec=22, char=SYN
{23},    // orig dec=23, char=ETB
{24},    // orig dec=24, char=CAN
{25},    // orig dec=25, char=EM
{26},    // orig dec=26, char=SUB
{27},    // orig dec=27, char=ESC
{28},    // orig dec=28, char=FS
{29},    // orig dec=29, char=GS
{30},    // orig dec=30, char=RS
{31},    // orig dec=31, char=US
{32},    // orig dec=32, char=(Space)
{33},    // orig dec=33, char=!
{34},    // orig dec=34, char="
{35},    // orig dec=35, char=#
{36},    // orig dec=36, char=$
{37},    // orig dec=37, char=%
{38},    // orig dec=38, char=&
{39},    // orig dec=39, char='
{40},    // orig dec=40, char=(
{41},    // orig dec=41, char=)
{42},    // orig dec=42, char=*
{43},    // orig dec=43, char=+
{44},    // orig dec=44, char=,
{45},    // orig dec=45, char=-
{46},    // orig dec=46, char=.
{47},    // orig dec=47, char=/
{48},    // orig dec=48, char=0
{49},    // orig dec=49, char=1
{50},    // orig dec=50, char=2
{51},    // orig dec=51, char=3
{52},    // orig dec=52, char=4
{53},    // orig dec=53, char=5
{54},    // orig dec=54, char=6
{55},    // orig dec=55, char=7
{56},    // orig dec=56, char=8
{57},    // orig dec=57, char=9
{58},    // orig dec=58, char=:
{59},    // orig dec=59, char=;
{60},    // orig dec=60, char=<
{61},    // orig dec=61, char==
{62},    // orig dec=62, char=>
{63},    // orig dec=63, char=?
{64},    // orig dec=64, char=@
{65},    // orig dec=65, char=A
{66},    // orig dec=66, char=B
{67},    // orig dec=67, char=C
{68},    // orig dec=68, char=D
{69},    // orig dec=69, char=E
{70},    // orig dec=70, char=F
{71},    // orig dec=71, char=G
{72},    // orig dec=72, char=H
{73},    // orig dec=73, char=I
{74},    // orig dec=74, char=J
{75},    // orig dec=75, char=K
{76},    // orig dec=76, char=L
{77},    // orig dec=77, char=M
{78},    // orig dec=78, char=N
{79},    // orig dec=79, char=O
{80},    // orig dec=80, char=P
{81},    // orig dec=81, char=Q
{82},    // orig dec=82, char=R
{83},    // orig dec=83, char=S
{84},    // orig dec=84, char=T
{85},    // orig dec=85, char=U
{86},    // orig dec=86, char=V
{87},    // orig dec=87, char=W
{88},    // orig dec=88, char=X
{89},    // orig dec=89, char=Y
{90},    // orig dec=90, char=Z
{91},    // orig dec=91, char=[
{92},    // orig dec=92, char=backslash
{93},    // orig dec=93, char=]
{94},    // orig dec=94, char=^
{95},    // orig dec=95, char=_
{96},    // orig dec=96, char=`
{65},    // orig dec=97, char=a
{66},    // orig dec=98, char=b
{67},    // orig dec=99, char=c
{68},    // orig dec=100, char=d
{69},    // orig dec=101, char=e
{70},    // orig dec=102, char=f
{71},    // orig dec=103, char=g
{72},    // orig dec=104, char=h
{73},    // orig dec=105, char=i
{74},    // orig dec=106, char=j
{75},    // orig dec=107, char=k
{76},    // orig dec=108, char=l
{77},    // orig dec=109, char=m
{78},    // orig dec=110, char=n
{79},    // orig dec=111, char=o
{80},    // orig dec=112, char=p
{81},    // orig dec=113, char=q
{82},    // orig dec=114, char=r
{83},    // orig dec=115, char=s
{84},    // orig dec=116, char=t
{85},    // orig dec=117, char=u
{86},    // orig dec=118, char=v
{87},    // orig dec=119, char=w
{88},    // orig dec=120, char=x
{89},    // orig dec=121, char=y
{90},    // orig dec=122, char=z
{123},    // orig dec=123, char={
{124},    // orig dec=124, char=|
{125},    // orig dec=125, char=}
{126},    // orig dec=126, char=~
{127},    // orig dec=127, char=DEL
{128},    // orig dec=128
{129},    // orig dec=129
{130},    // orig dec=130
{131},    // orig dec=131
{132},    // orig dec=132
{133},    // orig dec=133
{134},    // orig dec=134
{135},    // orig dec=135
{136},    // orig dec=136
{137},    // orig dec=137
{138},    // orig dec=138
{139},    // orig dec=139
{140},    // orig dec=140
{141},    // orig dec=141
{142},    // orig dec=142
{143},    // orig dec=143
{144},    // orig dec=144
{145},    // orig dec=145
{146},    // orig dec=146
{147},    // orig dec=147
{148},    // orig dec=148
{149},    // orig dec=149
{150},    // orig dec=150
{151},    // orig dec=151
{152},    // orig dec=152
{153},    // orig dec=153
{154},    // orig dec=154
{155},    // orig dec=155
{156},    // orig dec=156
{157},    // orig dec=157
{158},    // orig dec=158
{159},    // orig dec=159
{160},    // orig dec=160
{161},    // orig dec=161
{162},    // orig dec=162
{163},    // orig dec=163
{164},    // orig dec=164
{165},    // orig dec=165
{166},    // orig dec=166
{167},    // orig dec=167
{168},    // orig dec=168
{169},    // orig dec=169
{170},    // orig dec=170
{171},    // orig dec=171
{172},    // orig dec=172
{173},    // orig dec=173
{174},    // orig dec=174
{175},    // orig dec=175
{176},    // orig dec=176
{177},    // orig dec=177
{178},    // orig dec=178
{179},    // orig dec=179
{180},    // orig dec=180
{181},    // orig dec=181
{182},    // orig dec=182
{183},    // orig dec=183
{184},    // orig dec=184
{185},    // orig dec=185
{186},    // orig dec=186
{187},    // orig dec=187
{188},    // orig dec=188
{189},    // orig dec=189
{190},    // orig dec=190
{191},    // orig dec=191
{192},    // orig dec=192
{193},    // orig dec=193
{194},    // orig dec=194
{195},    // orig dec=195
{196},    // orig dec=196
{197},    // orig dec=197
{198},    // orig dec=198
{199},    // orig dec=199
{200},    // orig dec=200
{201},    // orig dec=201
{202},    // orig dec=202
{203},    // orig dec=203
{204},    // orig dec=204
{205},    // orig dec=205
{206},    // orig dec=206
{207},    // orig dec=207
{208},    // orig dec=208
{209},    // orig dec=209
{210},    // orig dec=210
{211},    // orig dec=211
{212},    // orig dec=212
{213},    // orig dec=213
{214},    // orig dec=214
{215},    // orig dec=215
{216},    // orig dec=216
{217},    // orig dec=217
{218},    // orig dec=218
{219},    // orig dec=219
{220},    // orig dec=220
{221},    // orig dec=221
{222},    // orig dec=222
{223},    // orig dec=223
{224},    // orig dec=224
{225},    // orig dec=225
{226},    // orig dec=226
{227},    // orig dec=227
{228},    // orig dec=228
{229},    // orig dec=229
{230},    // orig dec=230
{231},    // orig dec=231
{232},    // orig dec=232
{233},    // orig dec=233
{234},    // orig dec=234
{235},    // orig dec=235
{236},    // orig dec=236
{237},    // orig dec=237
{238},    // orig dec=238
{239},    // orig dec=239
{240},    // orig dec=240
{241},    // orig dec=241
{242},    // orig dec=242
{243},    // orig dec=243
{244},    // orig dec=244
{245},    // orig dec=245
{246},    // orig dec=246
{247},    // orig dec=247
{248},    // orig dec=248
{249},    // orig dec=249
{250},    // orig dec=250
{251},    // orig dec=251
{252},    // orig dec=252
{253},    // orig dec=253
{254},    // orig dec=254
{255},    // orig dec=255
};

