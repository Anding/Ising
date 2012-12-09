\ Ising model simulation in FORTH
\ Andrew Read
\ (C) 2012, GNU General Public License
\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\ functionality for running simulations with numerical results

8 constant n1						\ define the sample size as 2^n1, e.g. 256
5 constant n2						\ define the number of runs as 2^n2, e.g. 32
1 n1 lshift constant samples
1 32 n1 - lshift constant step
1 n2 lshift constant runcount
create magnetization{ samples 4 * allot
create energy{ samples 4 * allot
create heatcapacity{ samples 4 * allot

: } 
	4 * + 
;

: simulate ( Monte-Carlo Ising simulation, hard-coded paramarters for simplicity here)
	samples 0 DO					\ zero the results arrays
		0 magnetization{ i } !
		0 energy{ i } !
		0 heatcapacity{ i } !
	LOOP
	runcount 0 DO					\ simulation runs
		64 newlattice
		samples 0 DO			
			i step * Boltzmann !	\ step size
			1000000 run			\ equlibriate
			stats 				( C E M)
			abs magnetization{ i } +!
			energy{ i } +!
			heatcapacity{ i } +!
		LOOP	
	LOOP
	
;

: results ( output results.  Further need to divide by dimension * dimension * runcount to obtain per site values)
	cr ." Beta" 9 emit ." Magnetizaton" 9 emit ." Energy" 9 emit ." Heat capacity" cr
	samples 0 DO				
		i step * 10 u.r 9 emit
		magnetization{ i } @ 
		4 .r 9 emit
		energy{ i } @
		5 .r 9 emit 
		heatcapacity{ i } @
		4 .r cr
	LOOP	
;