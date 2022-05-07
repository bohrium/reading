/* author: samtenka
 * create: 2018-03-29 
 * change: 2018-03-29 
 * descrp: interfaces for ASCII Art translation, building on Bitmap 
 * usage : call `translate` and `stretched` 
*/

#ifndef ASCII_H
#define ASCII_H

#include "math.h"

#include "Bitmap.h"
#include <fstream>
#include <algorithm>
#include <iostream>

//const int types = 9;
//char asciis[] = " .!itmITM"; // 9 types
//int thresholds[] = {0, 32, 64, 96, 128, 160, 192, 224, 256};

const int types = 95;
char asciis[] = " !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"; 
int thresholds[95];

void init_thresholds(char const* in_name, int spacing) {
    Bitmap alpha;
    alpha.read_from(in_name);
    for (int i=0; i!=types; ++i) {
        thresholds[i]=0; 
    } 
    for (int c=0; c!=alpha.dims.width; ++c) {
        int i = c/spacing;
        for (int r=0; r!=alpha.dims.height; ++r) {
            thresholds[i] += 255 - alpha.data[r][c].R;
        }
    } 

    std::sort(asciis, asciis+types, [](int a, int b) {
        return thresholds[a-' '] < thresholds[b-' '];
    });

    int max = 0;
    for (int i=0; i!=types; ++i) {
        max = max < thresholds[i] ? thresholds[i] : max; 
    } 
    for (int i=0; i!=types; ++i) {
        thresholds[i] = (255 * thresholds[i])/max; 
    } 

    for (int i=0; i!=types; ++i) {
        std::cout << asciis[i] << "\t";
    } std::cout << std::endl; 
    for (int i=0; i!=types; ++i) {
        std::cout << thresholds[asciis[i]-' '] << "\t";
    } 

}

static int index(double value) {
    for(int i=0; i<types; ++i) {
        if(thresholds[asciis[i]-' ']>value) {
            return i;
        }
    }
    return types-1;
}

static double value(RGB rgb) {
    return 255 - (rgb.R + rgb.G + rgb.B)/3;
}

void translate(const Bitmap& bmp, char const* out_name ) {
    std::FILE* ascii = std::fopen(out_name, "w");

    double val_err = 0.0;
    for (int r=bmp.dims.height-1; r>=0; --r) {
        for (int c=0; c<bmp.dims.width; ++c) {
            val_err += value(bmp.data[r][c]);
            int ind = index(val_err);
            std::fputc(asciis[ind], ascii);
            val_err -= thresholds[asciis[ind]-' '];
        }
        std::fputc('\n', ascii);
    }

    fclose(ascii);
}

void stretch(const Bitmap& bmp, Bitmap& bmp2) {
    bmp2.allocate({bmp.dims.height/2, bmp.dims.width});
    for (int r=0; r!=bmp.dims.height/2; ++r) {
        for (int c=0; c!=bmp.dims.width; ++c) {
            bmp2.data[r][c].R = (bmp.data[2*r][c].R + bmp.data[2*r+1][c].R)/2;
            bmp2.data[r][c].G = (bmp.data[2*r][c].G + bmp.data[2*r+1][c].G)/2;
            bmp2.data[r][c].B = (bmp.data[2*r][c].B + bmp.data[2*r+1][c].B)/2;
        }
    }
}

typedef struct {
    float mean, variance;
} Stats;

Stats stats_from(const Bitmap& bmp, int r, int c, int d, float w) {
    float denom = 0.01;
          denom = 1       * denom;
    float sum   = 128     * denom;
    float sum2  = 128*128 * denom;
    for (int dr=-d; dr!=d+1; ++dr) {
        if (!(0<=r+dr && r+dr<bmp.dims.height)) { continue; }
        for (int dc=-d; dc!=d+1; ++dc) {
            if (!(0<=c+dc && c+dc<bmp.dims.width)) { continue; }
            float weight = 1.0/(w*w + dr*dr+dc*dc);
            float val = ((float)bmp.data[r+dr][c+dc].R +
                         (float)bmp.data[r+dr][c+dc].G +
                         (float)bmp.data[r+dr][c+dc].B  )/3.0; 
            denom += 1       * weight;
            sum   += val     * weight;
            sum2  += val*val * weight;
        }
    }
    sum /= denom;
    sum2 /= denom;
    return {sum, sum2-sum*sum};
} 

void contrast(const Bitmap& bmp, Bitmap& bmp2, float a, float b) {
    /* clear background from text and darken text 
     * 
     */
    bmp2.allocate({bmp.dims.height, bmp.dims.width});
    for (int r=0; r!=bmp.dims.height; ++r) {
        for (int c=0; c!=bmp.dims.width; ++c) {
            //Stats nine  = stats_from(bmp, r, c, 5, 3.0); 
            Stats seven = stats_from(bmp, r, c, 4, 2.4); 
            //Stats five  = stats_from(bmp, r, c, 3, 1.8); 
            //Stats three = stats_from(bmp, r, c, 2, 1.2); 
            Stats one   = stats_from(bmp, r, c, 1, 0.6); 

            if (seven.variance < 0.01*128*128) {
                bmp2.data[r][c].R = 255-24 + b*128*a;
                bmp2.data[r][c].G = 255-24 + b*128*a;
                bmp2.data[r][c].B = 255-24 + b*128*a;
            } else if (one.mean < seven.mean + a*sqrt(seven.variance)) {
                bmp2.data[r][c].R = 255-24 + b*128*(one.mean-seven.mean)/sqrt(seven.variance);
                bmp2.data[r][c].G = 255-24 + b*128*(one.mean-seven.mean)/sqrt(seven.variance);
                bmp2.data[r][c].B = 255-24 + b*128*(one.mean-seven.mean)/sqrt(seven.variance);
            } else {
                bmp2.data[r][c].R = 255-24 + b*128*a;
                bmp2.data[r][c].G = 255-24 + b*128*a;
                bmp2.data[r][c].B = 255-24 + b*128*a;
            }

#define CLIP(X) ((X)<0?0:(X)>255?255:(X))
            bmp2.data[r][c].R = CLIP(bmp2.data[r][c].R);
            bmp2.data[r][c].G = CLIP(bmp2.data[r][c].G);
            bmp2.data[r][c].B = CLIP(bmp2.data[r][c].B);
#define ROUND(X) ((X)<24?(X)*(X)/24.0:(X)>255-24?255-24+(255-(X))*(255-(X))/24.0:(X))
            bmp2.data[r][c].R = ROUND(bmp2.data[r][c].R);
            bmp2.data[r][c].G = ROUND(bmp2.data[r][c].G);
            bmp2.data[r][c].B = ROUND(bmp2.data[r][c].B);


        }
    }
}

#endif // ASCII_H

