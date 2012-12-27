\ Ising model simulation in FORTH
\ Andrew Read
\ (C) 2012, GNU General Public License
\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\ Graphical extensions for the NIGE Machine

: render-n ( display the lattice, for NIGE Machine) 
	csr-off 0 csr-x ! 0 csr-y !
	lattice
	dimension @ ROWS C@ 2 - min 0 DO
		dimension @ COLS c@ min 0 DO
			dup c@ 2 = if 159 else 160 then vemitraw
			1+
		LOOP 
		dimension @ COLS C@ - dup 0> IF 
			+			\ update for skipped columns
		ELSE	
			CR			\ newline
			drop
		THEN
	LOOP
	drop
	stats 
	cr
	." Magnetization " 4 .r 9 emit
	drop
	." Lattice energy " 5 .r 9 emit
	." Boltzmann " Boltzmann @ u.
;

: UMAX ( a b -- u, unsigned maximum of a and b)
	over over
	U> IF DROP ELSE NIP THEN
;

: UMIN ( a b -- u, unsigned minimum of a and b)
	over over
	U< IF DROP ELSE NIP THEN
;

: run-n ( -- indefinate iteration of the metropolis algorithm for NIGE Machine, press UP / DOWN to raise or lower temp, or ESC to stop)
	0 interlace cls				\ need non-interlaced screen
	BEGIN
		1000 0 DO metropolis LOOP
		render-n		
		key? IF
			key CASE
				4 OF Boltzmann dup @ 10000000 UMAX 1000000 - swap ! ENDOF
				5 OF Boltzmann dup @ 4284967295 UMIN 1000000 + swap ! ENDOF
				27 OF EXIT ENDOF
			ENDCASE
		THEN
	AGAIN
;
