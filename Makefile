
CFLAGS = -std=c99 -W -Wall -O3 -g
LDFLAGS = -lm

panorama: panorama.c

test: panorama
	./panorama 512x256 input.ppm

