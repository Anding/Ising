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
0 ,			\ expressed as x/8, in fractional format		
17044952 ,		\ a binary search algorithm finds bucket appropriate the given value of x/8 in this table
34648854 ,
52849611 ,
71689120 ,
91213850 ,
111475535 ,
132531995 ,
154448136 ,
177297155 ,
201162013 ,
226137257 ,
252331277 , 
279869154 ,
308896273 ,
339582962 ,
372130558 ,
406779413 ,
443819679 ,
483606094 ,
526578695 ,
573292572 ,
624461836 ,
681026832 ,
744261117 ,
815950238 ,
898709254 ,
996592395 ,
1116391676 ,
1270839813 ,
1488522235 ,
1860652794 ,
4294967295 ,


create exptab		\ evenly spaced probabilites in the range 0 <= x < 1
4294967295 ,		\ for each bucket in the table above, this table returns the mid-bucket probability given by Exp(-x)
4227858432 ,
4093640704 ,
3959422976 ,
3825205248 ,
3690987520 ,
3556769792 ,
3422552064 ,
3288334336 ,
3154116608 ,
3019898880 ,
2885681152 ,
2751463424 ,
2617245696 ,
2483027968 ,
2348810240 ,
2214592512 ,
2080374784 ,
1946157056 ,
1811939328 ,
1677721600 ,
1543503872 ,
1409286144 ,
1275068416 ,
1140850688 ,
1006632960 ,
872415232 ,
738197504 ,
603979776 ,
469762048 ,
335544320 ,
201326592 ,
67108864 ,

hex 
20000000 constant 1<<29 
80000000 constant 1<<31
decimal

: negexp ( frac int -- frac, return Exp[-x] where x is a postive number expressed as frac and int and int < 8)
		\ first convert frac int to a single frac by dividing by 8 and combining the fractional and integer parts
		1<<29 *				\ 29 lshift optimized using a multiply
		swap 
		3 rshift 
		or					\ x <- x/8, in fractional format (assumes 32 bit wordsize)
		\ use a binary search algorthm to find the appropriate bucket for x/8 in the log tab
		>R
		32 logtab 64 + 	( delta address)
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
1846835896 constant critical			\ critical beta for a 64*64 lattice, by experiment
create lattice 16384 allot				\ 128 * 128 lattice space allocation

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

: stats ( -- C E M, calculate the magnetization, total energy, and heat capacity of the lattice)
	0 0 0 					( E E2 M)
	modulus @  0 DO 
		i lattice + c@ 1- +		( E E2 M')
		rot				( E2 M' E)
		i energy >R R@ +		( E2 M E' R:energy loop)
		rot 				( M' E' E2 R:energy loop)
		R> dup * + 			( M' E' E2' R:loop)
		rot				( E' E2' M')
	LOOP 
	>R					( E E2 R:M)
	over dup * modulus @ /		( E E2 E^2 R:M)
	- 
	Boltzmann @ UM* nip		
	Boltzmann @ UM* nip			( E C R:M)
	swap R>				( C E M)
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

: simulate ( Monte-Carlo Ising simulation, hard-coded paramarters for simplicity here)
	cr ." Magnetization	Energy		Heat Capacity" cr
	64 newlattice
	4294967295 0 DO			\ range of Boltzmann beta values from 0.00000 to .99999
		i dup u. Boltzmann !
		1000000 run			\ equlibriate
		stats 9 emit . 9 emit . 9 emit  . cr
	42949672 +LOOP			\ step size through Boltzmann beta
;