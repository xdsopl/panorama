
CFLAGS = -std=c99 -W -Wall -O3 -D_GNU_SOURCE=1 -g
LDFLAGS = -lm

panorama: panorama.c

test: panorama
	./panorama 512x256 input.ppm
	convert output.ppm /var/www/localhost/htdocs/panorama.jpg

