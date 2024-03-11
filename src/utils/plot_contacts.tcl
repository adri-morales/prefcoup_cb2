set myfile "/home/adri/Documents/projects/gpcrdd/bias/results/contact_coefficients.txt"
set a [open $myfile]
set lines [split [read $a] "\n"]
close $a;                          # Saves a few bytes :-)
foreach line $lines {
    # do something with each line...
    lassign $line res1 res2 coef
    set sel [atomselect top "resid ${res1} and name CA"]
    lassign [$sel get {x y z}] coord1
    
    set sel [atomselect top "resid ${res2} and name CA"]
    lassign [$sel get {x y z}] coord2
    
    if {$coef > 0} {
	   draw color red
	} else {
	   draw color blue
	}

    draw cylinder $coord1 $coord2 radius [expr abs($coef)]
    
}