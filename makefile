make:
	gcc FFT.c -lm -o FFT -std=gnu99 -Wall -pedantic
	gcc CooleyTukey.c -lm -o CT -std=gnu99 -Wall -pedantic
	gcc DTFT.c -lm -o DTFT -std=gnu99 -Wall -pedantic
	gcc goodthomas.c -lm -o GT -std=gnu99 -Wall -pedantic 
