#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include "display.h"

double display(char endr,double dx,double dt) {
	
		char gnuplotparam[80];
		sprintf(gnuplotparam, "gnuplot  -e dt=%f -e dx=%f -e end=%s -p  \"script.txt\" ", dt, dx,endr);
		system(gnuplotparam);
}