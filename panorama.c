/*
panorama - downsample spherical panorama images
Written in 2014 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

float srgb(float v)
{
	float K0 = 0.03928f;
	float a = 0.055f;
	float phi = 12.92f;
	float gamma = 2.4f;
	return v <= K0 / phi ? v * phi : (1.0f + a) * powf(v, 1.0f / gamma) - a;
}

float linear(float v)
{
	float K0 = 0.03928f;
	float a = 0.055f;
	float phi = 12.92f;
	float gamma = 2.4f;
	return v <= K0 ? v / phi : powf((v + a) / (1.0f + a), gamma);
}

struct rgb {
	float r, g, b;
};

struct image {
	struct rgb *buffer;
	int width, height, total;
	char *name;
};

void delete_image(struct image *image)
{
	free(image->buffer);
	free(image);
}

struct image *new_image(char *name, int width, int height)
{
	struct image *image = malloc(sizeof(struct image));
	image->height = height;
	image->width = width;
	image->total = width * height;
	image->name = name;
	image->buffer = malloc(sizeof(struct rgb) * width * height);
	return image;
}

void downsample(struct image *output, struct image *input)
{
	int ow = output->width;
	int oh = output->height;
	int iw = input->width;
	int ih = input->height;
	struct rgb *ob = output->buffer;
	struct rgb *ib = input->buffer;
	for (int oj = 0; oj < oh; oj++) {
		int ij = (ih * oj) / oh;
		for (int oi = 0; oi < ow; oi++) {
			int ii = (iw * oi) / ow;
			ob[ow * oj + oi] = ib[iw * ij + ii];
		}
	}
}

struct image *read_ppm(char *name)
{
	FILE *file = fopen(name, "r");
	if (!file) {
		fprintf(stderr, "could not open \"%s\" file to read.\n", name);
		return 0;
	}
	if ('P' != fgetc(file) || '6' != fgetc(file)) {
		fprintf(stderr, "file \"%s\" not P6 image.\n", name);
		fclose(file);
		return 0;
	}
	int integer[3];
	struct image *image = 0;
	int c = fgetc(file);
	if (EOF == c)
		goto eof;
	for (int i = 0; i < 3; i++) {
		while ('#' == (c = fgetc(file)))
			while ('\n' != (c = fgetc(file)))
				if (EOF == c)
					goto eof;
		while ((c < '0') || ('9' < c))
			if (EOF == (c = fgetc(file)))
				goto eof;
		char str[16];
		for (int n = 0; n < 16; n++) {
			if (('0' <= c) && (c <= '9') && n < 15) {
				str[n] = c;
				if (EOF == (c = fgetc(file)))
					goto eof;
			} else {
				str[n] = 0;
				break;
			}
		}
		integer[i] = atoi(str);
	}
	if (!(integer[0]|integer[1]|integer[2])) {
		fprintf(stderr, "could not read image file \"%s\".\n", name);
		fclose(file);
		return 0;
	}
	if (integer[2] != 255) {
		fprintf(stderr, "cant read \"%s\", only 8 bit per channel SRGB supported at the moment.\n", name);
		fclose(file);
		return 0;
	}
	image = new_image(name, integer[0], integer[1]);
	for (int i = 0; i < image->total; i++) {
		int r = fgetc(file);
		int g = fgetc(file);
		int b = fgetc(file);
		if (EOF == r || EOF == g || EOF == b)
			goto eof;
		image->buffer[i] = (struct rgb){ linear(r / 255.0f), linear(g / 255.0f), linear(b / 255.0f) };
	}
	fclose(file);
	return image;
eof:
	fprintf(stderr, "EOF while reading from \"%s\".\n", name);
	fclose(file);
	delete_image(image);
	return 0;
}

int write_ppm(struct image *image)
{
	FILE *file = fopen(image->name, "w");
	if (!file) {
		fprintf(stderr, "could not open \"%s\" file to write.\n", image->name);
		return 0;
	}
	if (!fprintf(file, "P6 %d %d 255\n", image->width, image->height)) {
		fprintf(stderr, "could not write to file \"%s\".\n", image->name);
		fclose(file);
		return 0;
	}
	for (int i = 0; i < image->total; i++) {
		if (EOF == fputc(255.0f * srgb(image->buffer[i].r), file))
			goto eof;
		if (EOF == fputc(255.0f * srgb(image->buffer[i].g), file))
			goto eof;
		if (EOF == fputc(255.0f * srgb(image->buffer[i].b), file))
			goto eof;
	}
	fclose(file);
	return 1;
eof:
	fprintf(stderr, "EOF while writing to \"%s\".\n", image->name);
	fclose(file);
	return 0;
}

int main(int argc, char **argv)
{
	if (argc != 3) {
		fprintf(stderr, "usage: %s wxh input.ppm\n", argv[0]);
		return 1;
	}
	char *width = argv[1];
	char *seperator = strchr(argv[1], 'x');
	*seperator = 0;
	char *height = seperator + 1;
	struct image *output = new_image("output.ppm", atoi(width), atoi(height));
	struct image *input = read_ppm(argv[2]);
	if (!input)
		return 1;
	if (input->width < output->width || input->height < output->height) {
		fprintf(stderr, "output %dx%d must be smaller or equal to input %dx%d\n", output->width, output->height, input->width, input->height);
		return 1;
	}
	downsample(output, input);
	if (!write_ppm(output))
		return 1;
	return 0;
}
