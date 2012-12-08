\ Ising model simulation in FORTH
\ Andrew Read
\ (C) 2012, GNU General Public License
\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\ Graphical extensions for the NIGE Machine

: render-n ( display the lattice, for NIGE Machine)
	csr-off 0 csr-x ! 3 csr-y !
	lattice
	dimension @ 0 DO
		9 emit
		dimension @ 0 DO
			dup c@ 2 = if 159 else 160 then emit
			1+
		LOOP CR
	LOOP
	drop
	stats 
	cr
	." Magnetization " . ."   " cr
	." Lattice energy " . ."   " cr
	." Lattice heat capacity " . ."   " 
;

: run-n ( -- indefinate iteration of the metropolis algorithm for NIGE Machine, press any key to stop)
	cls
	BEGIN
		1000 0 DO metropolis LOOP
		render-n		
		key?
	UNTIL
	key drop cr
;