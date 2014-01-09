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

struct uv {
	float u, v;
};

struct xyz {
	float x, y, z;
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

struct uv uv_sphere(struct xyz v)
{
	return (struct uv) {
		0.5f + atan2f(v.z, v.x) / (2.0f * M_PI),
		acosf(v.y) / M_PI
	};
}

struct xyz xyz_sphere(struct uv v)
{
	return (struct xyz) {
		sinf(v.v * M_PI) * cosf((v.u - 0.5f) * 2.0f * M_PI),
		cosf(v.v * M_PI),
		sinf(v.v * M_PI) * sinf((v.u - 0.5f) * 2.0f * M_PI)
	};
}

struct xyz xyz_smul(float a, struct xyz v)
{
	return (struct xyz) { a * v.x, a * v.y, a * v.z };
}

struct rgb rgb_smul(float a, struct rgb v)
{
	return (struct rgb) { a * v.r, a * v.g, a * v.b };
}

struct uv uv_smul(float a, struct uv v)
{
	return (struct uv) { a * v.u, a * v.v };
}

struct rgb rgb_add(struct rgb a, struct rgb b)
{
	return (struct rgb) { a.r + b.r, a.g + b.g, a.b + b.b };
}

struct xyz xyz_add(struct xyz a, struct xyz b)
{
	return (struct xyz) { a.x + b.x, a.y + b.y, a.z + b.z };
}

float xyz_length(struct xyz v)
{
	return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

float uv_length(struct uv v)
{
	return sqrtf(v.u * v.u + v.v * v.v);
}

struct xyz xyz_normalize(struct xyz v)
{
	return xyz_smul(1.0f / xyz_length(v), v);
}

struct xyz xyz_orthogonal(struct xyz v)
{
	struct xyz ox = { 0.0f, -v.z, v.y };
	struct xyz oy = { v.z, 0.0f, -v.x };
	struct xyz oz = { -v.y, v.x, 0.0f };
	struct xyz o = fabsf(v.x) < fabsf(v.y) ? (fabsf(v.x) < fabsf(v.z) ? ox : oz) : (fabsf(v.y) < fabsf(v.z) ? oy : oz);
	return xyz_normalize(o);
}

struct xyz xyz_cross(struct xyz a, struct xyz b)
{
	return (struct xyz) {
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x
	};
}

float gauss(float x, float y, float radius)
{
	float sigma = radius / 3.0f;
	return radius ? expf(- (x * x + y * y) / (2.0f * sigma * sigma)) / (2.0f * M_PI * sigma * sigma) : 1.0f;
}

void downsample(struct image *output, struct image *input)
{
	int ow = output->width;
	int oh = output->height;
	int iw = input->width;
	int ih = input->height;
	int radius = fmaxf(iw / ow, ih / oh) / 2.0f;
	float delta = 1.0f / fmaxf(iw / 2.0f, ih);
	// fprintf(stderr, "%d %f\n", radius, delta);
	struct rgb *ob = output->buffer;
	struct rgb *ib = input->buffer;
	for (int oj = 0; oj < oh; oj++) {
		struct uv pol;
		pol.v = oj / (float)oh;
		int weight = fminf(8.0f, 1.0f / sinf(pol.v * M_PI));
		//fprintf(stderr, "row: % 4d weight: % 2d\n", oj, weight);
		for (int oi = 0; oi < ow; oi++) {
			pol.u = oi / (float)ow;
			struct xyz car = xyz_sphere(pol);
			struct xyz orth0 = xyz_orthogonal(car);
			struct xyz orth1 = xyz_cross(orth0, car);
			float kernel_sum = 0.0f;
			struct rgb rgb_sum = { 0.0f, 0.0f, 0.0f };
			for (int aj = -radius * weight; aj <= radius * weight; aj++) {
				for (int ai = -radius * weight; ai <= radius * weight; ai++) {
					struct xyz ac = xyz_add(xyz_smul(delta * ai, orth0), xyz_smul(delta * aj, orth1));
					struct uv ap = uv_sphere(xyz_normalize(xyz_add(car, ac)));
					int ii = iw * ap.u;
					int ij = ih * ap.v;
					float kernel = gauss(ai, aj, radius * weight);
					rgb_sum = rgb_add(rgb_sum, rgb_smul(kernel, ib[iw * ij + ii]));
					kernel_sum += kernel;
				}
			}
			ob[ow * oj + oi] = rgb_smul(1.0f / kernel_sum, rgb_sum);
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

