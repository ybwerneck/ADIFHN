#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include "display.h"

void saveGif(char* endf, char* endr, double dx, double dt) {
	char gnuplotparam[200];
	sprintf(gnuplotparam, "gnuplot  -e dt=%f -e dx=%f -e endf=\'%s\' -e endr=\'%s\' -p  \"scripts\\salvarGif.txt\" ", dt, dx, endf, endr);
	system(gnuplotparam);
}

void exibirGif(char* endf, double dx, double dt) {
	char gnuplotparam[200];
	sprintf(gnuplotparam, "gnuplot  -e dt=%f -e dx=%f -e end=\'%s\' -p  \"scripts\\exibirGif.txt\" ", dt, dx, endf);
	system(gnuplotparam);
}

void saveFoto(char* endf, char* endr, double dx, double dt) {
	char gnuplotparam[200];
	sprintf(gnuplotparam, "gnuplot  -e dt=%f -e dx=%f -e endf=\'%s\' -e endr=\'%s\' -p  \"scripts\\salvarFoto.txt\" ", dt, dx, endf, endr);
	system(gnuplotparam);

}
void exibirFoto(char* endf, double dx, double dt, double dtTarget) {
	char gnuplotparam[200];
	sprintf(gnuplotparam, "gnuplot  -e dt=%f -e dx=%f -e end=\'%s\' -e dtT=\'%d\'  -p  \"scripts\\exibirFoto.txt\" ", dt, dx, endf, dtTarget);
	system(gnuplotparam);
}