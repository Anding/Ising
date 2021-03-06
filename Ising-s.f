\ Ising model simulation in FORTH
\ Andrew Read
\ (C) 2012, GNU General Public License
\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\ functionality for running simulations with numerical results

1000 constant runcount				\ define the number of runs for averaging
8 constant n1						\ define the numer of steps across temperature as 2^n1, e.g. 256
1 n1 lshift constant samples
1 32 n1 - lshift constant step
create magnetization{ samples 4 * allot
create energy{ samples 4 * allot
create energy2{ samples 4 * allot

: } 
	4 * + 
;

: simulate ( Monte-Carlo Ising simulation, hard-coded paramarters for simplicity here)
	samples 0 DO					\ zero the results arrays
		0 magnetization{ i } !
		0 energy2{ i } !		
		0 energy{ i } !
	LOOP
	64 newlattice
	samples 0 DO			
		i step * Boltzmann !		\ step size
		100000 run			\ pre-equlibriate
		runcount 0 DO
			50000 run
			stats 				( C E M)
			abs magnetization{ j } +!
			energy2{ j } +!			
			energy{ j } +!
		LOOP
	LOOP	
	
;

: results ( output results.  Further need to divide by dimension * dimension * runcount to obtain per site values)
	cr ." Beta" 9 emit ." Magnetizaton" 9 emit ." Energy squared" 9 emit ." Energy" cr
	samples 0 DO				
		i step * 10 u.r 9 emit
		magnetization{ i } @ 
		4 .r 9 emit
		energy2{ i } @
		5 .r 9 emit		
		energy{ i } @
		5 .r cr

	LOOP	
;