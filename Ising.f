variable Xn

: seed ( n --, seed the random series with n)
	Xn !
;

: rnd ( -- n, random number between 0 and 2^32 - 1)
	Xn dup @ 3141592621 UM*				\ Algorithm from Knuth 3.6
	1 0 D+
	drop
	dup rot !
;

create logtab		\ Logarithms corresponding to the exponentials in the table below
34648854 ,		\ expressed as x/8, in fractional format
71689120 ,
111475535 ,
154448136 ,
201162013 ,
252331277 ,
308896273 ,
372130558 ,
443819679 ,
526578695 ,
624461836 ,
744261117 ,
898709254 ,
1116391676 ,
1488522235 ,
4294967295 ,

create exptab		\ sixeteen evenly spaced values of x in the range 0 <= x < 1
4160749568 ,		\ expressed in fractional format
3892314112 ,
3623878656 ,
3355443200 ,
3087007744 ,
2818572288 ,
2550136832 ,
2281701376 ,
2013265920 ,
1744830464 ,
1476395008 ,
1207959552 ,
939524096 ,
671088640 ,		\ 5/32
402653184 ,		\ 3/32
134217728 ,		\ 1/32

hex 
08000000 constant 1<<27 
80000000 constant 1<<31
decimal

: negexp ( frac int -- frac, return Exp[-x] where x is a postive number expressed as frac and int)
		27 lshift		\ optimize with multiply
		swap 
		3 rshift 
		or			\ x <- x/8, in fractional format (assumes 32 bit wordsize)
		dup u. >R
		0 BEGIN					\ find the lowest log which x exceeds
			dup logtab + @
			R@ U> not
		WHILE
			4 +
		REPEAT
		R> drop
		exptab + @					\ lookup the corresponding exponential
;
	
	
variable dimension
variable modulus
variable lattice

: newlattice ( n - prepare an n * n lattice)
	dup dimension !			( n)
	dup * dup modulus !			( nn)
	dup allocate				( nn addr ior) \ obtain memory
	if
		." out of memory"
		drop exit
	then
	dup lattice !				( nn addr)
	dup rot + swap			( nn end-addr addr)
	DO							\ randomize the lattice
		rnd 1<<31 and 				\ msb
		0= i c!					\ write byte value
	LOOP
; 

: render ( display the lattice)
	cr
	lattice @
	dimension @ 0 DO
		dimension @ 0 DO
			dup c@ if ." + " else ." - " then
			1+
		LOOP CR
	LOOP
	drop
;

: rnd-cell ( -- n, pick a random cell in the lattice)
	modulus @ rnd um* nip
;

: neighbours ( n -- n1 n2 n3 n4, find the nearest neighbours of the cell)
	modulus @ 1- >R 			\ circular boundary conditions
	1+ R@ and					( right)
	dup 2 - R@ and				( right left)
	dup 1+ dimension @ dup rot + R@ and	( right left dim down)
	dup rot 2* - R> and				( right left down up)
;
	