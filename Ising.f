\ Ising model simulation in FORTH
\ Andrew Read
\ (C) 2012, GNU General Public License
\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\ Random number generator - algorithm from Kunth 3.6
variable Xn		\ last random number

: seed ( u --, seed the random series with n)
	Xn !
;

: rnd ( -- u, random number between 0 and 2^32 - 1, often treated as a fraction)
	Xn dup @ 3141592621 UM*				\  3141592621
	1 0 D+
	drop
	dup rot !
;

: t rnd 32 u.r ;

: rnd-int ( u -- u, pick a random integer between 0 and u)
	rnd um* nip
;

\ Special purpose lookup tables for calculating Exp( - Beta * DeltaE)

create logtab		\ Logarithms corresponding to the exponentials in the table below	
			\ expressed as x/8, in fractional format	
			\ a binary search algorithm finds bucket appropriate the given value of x/8 in this table
0 ,
67108864 ,
100663296 ,
134217728 ,
167772160 ,
201326592 ,
234881024 ,
268435456 ,
301989888 ,
335544320 ,
369098752 ,
402653184 ,
436207616 ,
469762048 ,
503316480 ,
536870912 ,
570425344 ,
603979776 ,
637534208 ,
671088640 ,
704643072 ,
738197504 ,
771751936 ,
805306368 ,
838860800 ,
872415232 ,
905969664 ,
939524096 ,
973078528 ,
1006632960 ,
1040187392 ,
1073741824 ,
1107296256 ,
1140850688 ,
1174405120 ,
1207959552 ,
1241513984 ,
1275068416 ,
1308622848 ,
1342177280 ,
1375731712 ,
1409286144 ,
1442840576 ,
1476395008 ,
1509949440 ,
1543503872 ,
1577058304 ,
1610612736 ,
1644167168 ,
1677721600 ,
1711276032 ,
1744830464 ,
1778384896 ,
1811939328 ,
1845493760 ,
1879048192 ,
1912602624 ,
1946157056 ,
1979711488 ,
2013265920 ,
2046820352 ,
2080374784 ,
2113929216 ,
2147483648 ,
2181038080 ,
2214592512 ,
2248146944 ,
2281701376 ,
2315255808 ,
2348810240 ,
2382364672 ,
2415919104 ,
2449473536 ,
2483027968 ,
2516582400 ,
2550136832 ,
2583691264 ,
2617245696 ,
2650800128 ,
2684354560 ,
2717908992 ,
2751463424 ,
2785017856 ,
2818572288 ,
2852126720 ,
2885681152 ,
2919235584 ,
2952790016 ,
2986344448 ,
3019898880 ,
3053453312 ,
3087007744 ,
3120562176 ,
3154116608 ,
3187671040 ,
3221225472 ,
3254779904 ,
3288334336 ,
3321888768 ,
3355443200 ,
3388997632 ,
3422552064 ,
3456106496 ,
3489660928 ,
3523215360 ,
3556769792 ,
3590324224 ,
3623878656 ,
3657433088 ,
3690987520 ,
3724541952 ,
3758096384 ,
3791650816 ,
3825205248 ,
3858759680 ,
3892314112 ,
3925868544 ,
3959422976 ,
3992977408 ,
4026531840 ,
4060086272 ,
4093640704 ,
4127195136 ,
4160749568 ,
4194304000 ,
4227858432 ,
4261412864 ,
4294967295 ,


create exptab		\ evenly spaced probabilites in the range 0 <= x < 1
			\ for each bucket in the table above, this table returns the mid-bucket probability given by Exp(-x)
4162825044 ,
3910612224 ,
3673680207 ,
3451103175 ,
3242011404 ,
3045587862 ,
2861065022 ,
2687721855 ,
2524881020 ,
2371906212 ,
2228199679 ,
2093199885 ,
1966379315 ,
1847242415 ,
1735323655 ,
1630185710 ,
1531417750 ,
1438633839 ,
1351471421 ,
1269589907 ,
1192669343 ,
1120409161 ,
1052527001 ,
988757614 ,
928851818 ,
872575531 ,
819708853 ,
770045204 ,
723390523 ,
679562507 ,
638389896 ,
599711808 ,
563377106 ,
529243813 ,
497178551 ,
467056025 ,
438758531 ,
412175496 ,
387203045 ,
363743598 ,
341705488 ,
321002599 ,
301554034 ,
283283799 ,
266120501 ,
249997075 ,
234850518 ,
220621644 ,
207254855 ,
194697918 ,
182901767 ,
171820309 ,
161410243 ,
151630891 ,
142444039 ,
133813791 ,
125706423 ,
118090256 ,
110935529 ,
104214285 ,
97900261 ,
91968784 ,
86396677 ,
81162167 ,
76244800 ,
71625361 ,
67285800 ,
63209159 ,
59379510 ,
55781887 ,
52402233 ,
49227342 ,
46244809 ,
43442977 ,
40810900 ,
38338293 ,
36015493 ,
33833425 ,
31783561 ,
29857892 ,
28048894 ,
26349497 ,
24753062 ,
23253350 ,
21844501 ,
20521009 ,
19277704 ,
18109727 ,
17012514 ,
15981778 ,
15013491 ,
14103869 ,
13249359 ,
12446621 ,
11692518 ,
10984104 ,
10318611 ,
9693438 ,
9106142 ,
8554429 ,
8036142 ,
7549257 ,
7091871 ,
6662196 ,
6258554 ,
5879367 ,
5523154 ,
5188523 ,
4874167 ,
4578856 ,
4301437 ,
4040826 ,
3796005 ,
3566016 ,
3349962 ,
3146998 ,
2956331 ,
2777216 ,
2608953 ,
2450885 ,
2302393 ,
2162898 ,
2031855 ,
1908751 ,
1793105 ,
1684467 ,
1582410 ,
1486536 ,

1 24 lshift constant 1<<24 
1 29 lshift constant 1<<29 
1 31 lshift constant 1<<31


: negexp ( frac int -- frac, return Exp[-x] where x is a postive number expressed as frac and int and int < 8)
		\ first convert frac int to a single frac by dividing by 8 and combining the fractional and integer parts
		1<<29 *				\ 29 lshift optimized using a multiply
		swap 
		3 rshift 
		or					\ x <- x/8, in fractional format (assumes 32 bit wordsize)
		\ use a binary search algorthm to find the appropriate bucket for x/8 in the log tab
		>R
		128 logtab 64 + 	( delta address)   \ delta = (tablesize * 4 bytes) / 4 quarters
		BEGIN				
			dup @
			R@ U< IF			\ bucket is smaller than x/8, continue search
				over +			\ add the delta to address to move up the table
				true			\ flag to indicate repeat
			ELSE				\ bucket is greater than x/8
				dup 4 - @
				R@ U> IF		\ prior bucket is also greater than than x/8, continue search
					over -
					true				
				else			\ prior bucket is smaller than or equal to x/8, so this is the correct bucket
					false		\ flag to indicate end search		
				THEN
			THEN
		WHILE
			swap 2/			\ halve the value of delta for the next iteration
			4 max				\ minimum value of delta must be 4 - the size of a table entry
			swap		( delta' address')
		REPEAT
		nip R> drop		( address)
		logtab - exptab + @			\ lookup the corresponding exponential
;

	
variable dimension					\ dimension is the width and height of the array
variable modulus					\ modulus = width*height
variable Boltzmann					\ Boltzmann beta ( 1 / T )
1892735628 constant critical			\ critical beta, 2^32 / 2.269185
create lattice 16384 allot				\ lattice space allocation sufficient for 128 * 128 

: newlattice ( n - prepare an n * n lattice, initiation of an Ising model experiment, BUT n MUST BE A POWER OF TWO)
	dup dimension !			( n)
	dup * dup modulus !			( nn)
	lattice dup rot + swap			( nn end-addr addr)
	DO							\ randomize the lattice
		rnd 1<<31 and IF				\ msb
			2						
		ELSE
			0
		THEN
		i c!						\ write byte value
	LOOP
	4294967295 Boltzmann !				\ reset temperature
; 

: rnd-cell ( -- n, pick a random cell in the lattice)
	modulus @ rnd-int
;

: neighbours ( n -- n1 n2 n3 n4, find the nearest neighbours of the cell)
	modulus @ 1- >R 			\ circular boundary conditions BUT dimension MUST BE A POWER OF TWO
	1+ R@ and					( right)
	dup 2 - R@ and				( right left)
	dup 1+ dimension @ dup rot + R@ and	( right left dim down)
	dup rot 2* - R> and				( right left down up)
;

: sum-neighbours ( n1 n2 n3 n4 -- m, sum the states of the neighbours)
	lattice >R
	R@ + C@
	swap R@ + C@ +
	swap R@ + C@ +
	swap R> + C@ +
	4 -
;

: energy ( n -- e, sum the interaction energy attributable to cell n)
	dup lattice + c@ 1- swap		( cell-m cell)
	neighbours sum-neighbours 		( cell-m neigh-m)
	* invert 1+				\ multiply and then negate
;

: stats ( -- E E2 M, calculate the magnetization,  sum of energy, and sum of energy squared)
	0 0 0 					( E E2 M)
	modulus @  0 DO 
		i lattice + c@ 1- +		( E E2 M')
		rot				( E2 M' E)
		i energy >R R@ +		( E2 M E' R:energy loop)
		rot 				( M' E' E2 R:energy loop)
		R> dup * + 			( M' E' E2' R:loop)
		rot				( E' E2' M')
	LOOP 
	\ >R					( E E2 R:M)
	\ over dup * modulus @ /		( E E2 E^2 R:M)
	\ - 
	\ Boltzmann @ UM* nip		
	\ Boltzmann @ UM* nip		( E C R:M)
	\ swap R>				( C E M)
;

: flip ( n --, flip cell n)
	lattice + dup				( addr addr)
	c@ 2 swap -				( addr flip)
	swap c!
;
	
: metropolis ( -- select and flip a cell according to the Metropolis rule)
	rnd-cell dup energy			( cell energy --)
	dup 0< IF					\ a flip will raise the energy of the system
		invert 1+ 2*			( cell delta-energy)
		Boltzmann @ UM* negexp	( cell p)
		rnd				( cell p rnd)
		U> 				( cell flag)
	ELSE						\ a flip will lower the energy of the system, or leave unchanged
		drop -1			( cell true)
	THEN
	IF flip ELSE drop THEN			\ conditional flip
;

: run ( n -- run n iterations of the metropolis algorithm, used to equlibriate at a given Boltzmann beta)
	0 DO metropolis LOOP				\ need to (1) initiate with newlattice and (2) set Boltzmann beta, first
;

: render ( display the lattice in ASCII format)
	cr						
	lattice
	dimension @ 0 DO
		dimension @ 0 DO
			dup c@ 2 = if ." + " else ." - " then
			1+
		LOOP CR
	LOOP
	drop
;