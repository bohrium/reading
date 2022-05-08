/* author: samtenka
 * create: 2022-05-07 
 * change: 2018-03-29 
 * descrp: interfaces for ASCII Art translation, building on Bitmap 
 * usage : call `translate` and `clear_background` 
*/

#ifndef ASCII_H
#define ASCII_H

#include "math.h"

#include "Bitmap.h"
#include <fstream>
#include <algorithm>
#include <iostream>

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
    int locs[4][2];
} FourDotLocations;

FourDotLocations find_reds(const Bitmap& bmp)
{
    /* find four red dots */

    FourDotLocations reds;
    int nb_reds_found=0;
    for (int r=0; r!=bmp.dims.height && nb_reds_found!=4; ++r) {
        for (int c=0; c!=bmp.dims.width && nb_reds_found!=4; ++c) {
            if (!(bmp.data[r][c].R == 255 &&   
                  bmp.data[r][c].G ==   0 &&  
                  bmp.data[r][c].B ==   0   )) { continue; }
            bool too_close = false;
            for (int n=0; n!=nb_reds_found; ++n) {
#define EUCLID2(R,RR,C,CC) (((R)-(RR))*((R)-(RR))+((C)-(CC))*((C)-(CC)))
                if (EUCLID2(r,reds.locs[n][0],c,reds.locs[n][1]) < 50*50) { too_close=true; }
            }
            if (too_close) { continue; }
            reds.locs[nb_reds_found][0] = r;
            reds.locs[nb_reds_found][1] = c;
            nb_reds_found += 1;
        }
    }
    if (nb_reds_found!=4) { std::cout << "\n:(\n" << std::flush; return reds; }

    /* sort */
#define SWAP(X,Y) { int tmp; tmp=(X); (X)=(Y); (Y)=tmp; }
#define ALIGN_ROW(i,j) {if (reds.locs[i][0]>reds.locs[j][0]) { SWAP(reds.locs[i][0], reds.locs[j][0]); SWAP(reds.locs[i][1], reds.locs[j][1]); }}
#define ALIGN_COL(i,j) {if (reds.locs[i][1]>reds.locs[j][1]) { SWAP(reds.locs[i][0], reds.locs[j][0]); SWAP(reds.locs[i][1], reds.locs[j][1]); }}

    ALIGN_ROW(0,1); ALIGN_ROW(1,2); ALIGN_ROW(2,3);
    ALIGN_ROW(0,1); ALIGN_ROW(1,2);
    ALIGN_ROW(0,1);
    ALIGN_COL(0,1);
    ALIGN_COL(2,3);

    std::cout << "\n " << reds.locs[0][0] << " : " << reds.locs[0][1] << "\n"; 
    std::cout << "\n " << reds.locs[1][0] << " : " << reds.locs[1][1] << "\n"; 
    std::cout << "\n " << reds.locs[2][0] << " : " << reds.locs[2][1] << "\n"; 
    std::cout << "\n " << reds.locs[3][0] << " : " << reds.locs[3][1] << "\n"; 

    return reds;
}

void reframe(const Bitmap& bmp, Bitmap& bmp2)
{
    FourDotLocations reds = find_reds(bmp);

    int height = (reds.locs[2][0]+reds.locs[3][0] - reds.locs[1][0]-reds.locs[0][0])/2; 
    int width  = (reds.locs[1][1]+reds.locs[3][1] - reds.locs[0][1]-reds.locs[2][1])/2; 
    bmp2.allocate({height, width});

    for (int r=0; r!=bmp2.dims.height; ++r) {
        for (int c=0; c!=bmp2.dims.width; ++c) {
            float y = ((float)r)/bmp2.dims.height; 
            float x = ((float)c)/bmp2.dims.width; 
            int rr = (1-y)*((1-x)*(float)reds.locs[0][0]+x*(float)reds.locs[1][0]) + y*((1-x)*(float)reds.locs[2][0]+x*(float)reds.locs[3][0]); 
            int cc = (1-y)*((1-x)*(float)reds.locs[0][1]+x*(float)reds.locs[1][1]) + y*((1-x)*(float)reds.locs[2][1]+x*(float)reds.locs[3][1]); 
            if ((0<=rr && rr<bmp.dims.height) && 
                (0<=cc && cc<bmp.dims.width )   ) {
                bmp2.data[r][c].R = bmp.data[rr][cc].R;
                bmp2.data[r][c].G = bmp.data[rr][cc].G;
                bmp2.data[r][c].B = bmp.data[rr][cc].B;
            }
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

float const BACKGROUND_MARGIN = 24.0;
float const ROUND_MARGIN = 16.0;
float const BACKGROUND_VARIANCE = 0.01;
void clear_background(const Bitmap& bmp, Bitmap& bmp2, float stdthresh, float darken_factor) {
    /* clear background from text and darken text 
     */
    bmp2.allocate({bmp.dims.height, bmp.dims.width});
    for (int r=0; r!=bmp.dims.height; ++r) {
        for (int c=0; c!=bmp.dims.width; ++c) {
            Stats wide = stats_from(bmp, r, c, 4, 2.4); 
            Stats narr  = stats_from(bmp, r, c, 1, 0.6); 

            if (wide.variance < BACKGROUND_VARIANCE*128*128) {
                bmp2.data[r][c].R = 255-BACKGROUND_MARGIN + darken_factor*128*stdthresh;
                bmp2.data[r][c].G = 255-BACKGROUND_MARGIN + darken_factor*128*stdthresh;
                bmp2.data[r][c].B = 255-BACKGROUND_MARGIN + darken_factor*128*stdthresh;
            } else {
                float rel_brightness = (narr.mean - wide.mean) / sqrt(wide.variance);
#define MIN(X,Y) ((X)<(Y)?(X):(Y))
                rel_brightness = MIN(+1, rel_brightness); 
                bmp2.data[r][c].R = 255-BACKGROUND_MARGIN + darken_factor*128*rel_brightness;
                bmp2.data[r][c].G = 255-BACKGROUND_MARGIN + darken_factor*128*rel_brightness;
                bmp2.data[r][c].B = 255-BACKGROUND_MARGIN + darken_factor*128*rel_brightness;
            }

#define CLIP(X) ((X)<0?0:(X)>255?255:(X))
            bmp2.data[r][c].R = CLIP(bmp2.data[r][c].R);
            bmp2.data[r][c].G = CLIP(bmp2.data[r][c].G);
            bmp2.data[r][c].B = CLIP(bmp2.data[r][c].B);
#define ROUND(X) ((X)<ROUND_MARGIN?(X)*(X)/ROUND_MARGIN:(X)>255-ROUND_MARGIN?255-ROUND_MARGIN+(255-(X))*(255-(X))/ROUND_MARGIN:(X))
            bmp2.data[r][c].R = ROUND(bmp2.data[r][c].R);
            bmp2.data[r][c].G = ROUND(bmp2.data[r][c].G);
            bmp2.data[r][c].B = ROUND(bmp2.data[r][c].B);

        }
    }
}

#endif // ASCII_H

