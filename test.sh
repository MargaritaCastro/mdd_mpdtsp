#!/bin/bash

instance="../class3/m10Q10s111.tsp"
restuls="test.txt"

# Example :
# -c3 		: Medium-size instance of class 3
# -w1024 	: MDD of width 1024
# -d   		: DepthFirstSearch = TRUE - recommended for MDD implementation
# -rc 		: Capacity-based refinement
# -lt 		: Lagrangian penalties based on tour constrains

./mdd_mPDTSP $instance -c3  -w1024 -d -lt -rc -o$results
