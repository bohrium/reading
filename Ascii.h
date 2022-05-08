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

int const nb_types = 95;
char asciis[nb_types+1] = " !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"; 

int const max_char_height = 48; 
int const max_char_width  = 32; 

int const nb_fonts =  7; 

int ascii_heights[nb_fonts]          ;
int ascii_widths [nb_fonts][nb_types];
float ascii_shapes[nb_fonts][nb_types][max_char_height][max_char_width]; // mostly 0.0s; sparse 1.0s (so inverted)

bool is_light_row(const Bitmap& bmp, int r)
{
    int darks = 0;
    for (int c=0; c!=bmp.dims.width; ++c) {
        if (bmp.data[r][c].R < 255-16) { darks += 1; }
    }
    return darks < bmp.dims.width/15.0;
}
int next_light_row(const Bitmap& bmp, int r_start)
{
    int r = r_start;
    while (r!=bmp.dims.height && !is_light_row(bmp, r)) { ++r; }
    return r;
}
int next_nonlight_row(const Bitmap& bmp, int r_start)
{
    int r = r_start;
    while (r!=bmp.dims.height && is_light_row(bmp, r)) { ++r; }
    return r;
}

bool is_white_row(const Bitmap& bmp, int r)
{
    int darks = 0;
    for (int c=0; c!=bmp.dims.width; ++c) {
        if (bmp.data[r][c].R < 255-16) { darks += 1; }
    }
    return darks < bmp.dims.width/15.0; /* ATTENTION! to 15.0 should be 35.0 */
}
int next_white_row(const Bitmap& bmp, int r_start)
{
    int r = r_start;
    while (r!=bmp.dims.height && !is_white_row(bmp, r)) { ++r; }
    return r;
}
int next_nonwhite_row(const Bitmap& bmp, int r_start)
{
    int r = r_start;
    while (r!=bmp.dims.height && is_white_row(bmp, r)) { ++r; }
    return r;
}

bool is_red_col(const Bitmap& bmp, int r_start, int r_end, int c)
{
    for (int r=r_start; r!=r_end; ++r) {
        if (bmp.data[r][c].R > bmp.data[r][c].G+32) { return true; }
    }
    return false;
}
int next_red_col(const Bitmap& bmp, int r_start, int r_end, int c_start)
{
    int c = c_start;
    while (c!=bmp.dims.width && !is_red_col(bmp, r_start, r_end, c)) { ++c; }
    return c;
}
int next_nonred_col(const Bitmap& bmp, int r_start, int r_end, int c_start)
{
    int c = c_start;
    while (c!=bmp.dims.width && is_red_col(bmp, r_start, r_end, c)) { ++c; }
    return c;
}

bool is_white_col(const Bitmap& bmp, int r_start, int r_end, int c)
{
    for (int r=r_start; r!=r_end; ++r) {
        if (bmp.data[r][c].R <192) { return false; }
    }
    return true;
}
int next_white_col(const Bitmap& bmp, int r_start, int r_end, int c_start)
{
    int c = c_start;
    while (c!=bmp.dims.width && !is_white_col(bmp, r_start, r_end, c)) { ++c; }
    return c;
}
int next_nonwhite_col(const Bitmap& bmp, int r_start, int r_end, int c_start)
{
    int c = c_start;
    while (c!=bmp.dims.width && is_white_col(bmp, r_start, r_end, c)) { ++c; }
    return c;
}



void load_fonts(char const* in_name) {
    Bitmap alpha;
    alpha.read_from(in_name);
    int r_start, r_end;
    int c_start, c_end;

    r_end = 0;
    for (int f=0; f!=nb_fonts; ++f) {
        r_start = next_nonwhite_row(alpha, r_end  );
        r_end   = next_white_row   (alpha, r_start);

        //c_end = 0;
        //for (int t=0; t!=nb_types; ++t) {
        //    c_end   = next_red_col     (alpha, r_start, r_end, c_end  );
        //    c_start = next_nonred_col  (alpha, r_start, r_end, c_end  );
        //    c_start = next_nonwhite_col(alpha, r_start, r_end, c_start);
        //    c_end   = next_white_col   (alpha, r_start, r_end, c_start);
        
        ascii_heights[f] = r_end-r_start;
        std::cout << " " << ascii_heights[f] << " : ";

        c_end = next_red_col     (alpha, r_start, r_end, 0);
        for (int t=0; t!=nb_types; ++t) {
            c_start = next_nonred_col  (alpha, r_start, r_end, c_end  );
            c_end   = next_red_col     (alpha, r_start, r_end, c_start);

            ascii_widths [f][t] = c_end-c_start;
            std::cout << ascii_widths [f][t] << ", ";
            for (int r=0; r!=r_end-r_start; ++r) {
                for (int c=0; c!=c_end-c_start; ++c) {
                    ascii_shapes[f][t][r][c] = (255 - alpha.data[r_start+r][c_start+c].G)/255.0; 
                }
            }
        }
        std::cout << std::endl;
        std::cout << std::endl;
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

float match(const Bitmap& bmp, int r, int c, int f, int t, float scale_r, float scale_c)
{
    float vw = 0.0; 
    float ww = 0.0;
    float vv = 0.0;
    //float score = 0.0; 
    //float denom = 0.0;
    for (int dr=0; dr!=ascii_heights[f]; ++dr) {
        for (int dc=0; dc!=ascii_widths [f][t]; ++dc) {
            float v =  0.01 + (255 - bmp.data[r+(int)(dr*scale_r)][c+(int)(dc*scale_c)].G)/255.0;
            float w =  0.01 + ascii_shapes[f][t][dr][dc];
            //score += ((v>0.6&&w>0.6) ? 1.0 : 0.0);
            //score += ((v<0.6&&w<0.6) ? 0.8 : 0.0);
            //denom += ((       w>0.6) ? 1.0 : 0.0);
            //denom += ((       w<0.6) ? 0.8 : 0.0);

            vw += v*w;
            ww += w*w;
            vv += v*v;
        }
    }
    //return denom ? score/denom : 0.0;
    return (vv*ww ? vw/sqrt(ww*vv) : 0.0); 
    //return (vv*ww ? vw/sqrt(ww*ascii_heights[f]*ascii_widths[f][t]) : 0.0); 
}

typedef struct {
    int best_key;
    int c_offset;
} MatchType; 

int const r_jitter = 2;//+2;
int const c_jitter = 8;//-2;
MatchType best_match(const Bitmap& bmp, int r, int c, float line_height)
{
    float best_val = -2.0;
    int   best_key = -1  ;
    int   best_c_offset = 0;
    for (int f =1 ; f!= 7; ++f) {
        for (int t = 1; t!=nb_types; ++t) {
        //for (int t = 'A'-' '; t!='z'-' '; ++t) {
            if ('['-' '<=t && t<='`'-' ') { continue; }
            if (':'-' '<=t && t<='@'-' ') { continue; }
            if ('!'-' '<=t && t<='/'-' ') { continue; }
            if ('{'-' '<=t && t<='~'-' ') { continue; }
            for (int dr=-r_jitter; dr<r_jitter+1; dr+=1) {
                for (int dc=0        ; dc<c_jitter+1; dc+=1) {
                    float sfs_r[3] = {0.95,1.00,1.05};
                    //float sfs_c[6] = {0.82,0.91,1.00,1.09,1.18,1.27};
                    //float sfs_c[6] = {0.70,0.85,1.00,1.15,1.30,1.45};
                    float sfs_c[5] = {0.85, 0.95,1.00,1.05, 1.15};
                    for (int sfidx_r=1; sfidx_r!=2; ++sfidx_r) {
                        for (int sfidx_c=0; sfidx_c!=5; ++sfidx_c) {
                            float scale = line_height / ascii_heights[f];
                            float m = match(bmp, r+dr, c+dc, f, t, sfs_r[sfidx_r]*scale, sfs_c[sfidx_c]*scale);
                            if ( m <= best_val ) { continue; } 
                            best_val = m; 
                            best_key = t;
                            best_c_offset = dc+ascii_widths[f][t]*sfs_c[sfidx_c]*scale;
                            best_c_offset = best_c_offset<=0 ? 1 : best_c_offset;  
                        }
                    }
                }
            }
        }
    }
    if (best_val > 0.9 ) {
        return {best_key, best_c_offset};
    } else {
        return {'.'-' ', best_c_offset};
    }
}

void ocr(const Bitmap& bmp)
{
    int r = 0;
    while (true) {
        r = next_nonlight_row(bmp, r);
        if (!(r+r_jitter+max_char_height<=bmp.dims.height)) { break; } 
        int line_height = next_white_row(bmp, r) - r;
        std::cout << ":" << r << ":" << r+line_height << std::endl;

        int c = next_nonwhite_col(bmp, r, r+line_height, 0);
        while (c+c_jitter+max_char_width<=bmp.dims.width) {
            //std::cout << " " << c << ", " << std::flush;
            MatchType mt = best_match(bmp, r, c, line_height); 
            std::cout << (char)(' '+mt.best_key) << std::flush; 
            c += mt.c_offset;
        }
        std::cout << std::endl;

        r = next_light_row(bmp, r);
        if (!(r+r_jitter+max_char_height<=bmp.dims.height)) { break; } 
    } 
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




//static int index(double value) {
//    for(int i=0; i<types; ++i) {
//        if(thresholds[asciis[i]-' ']>value) {
//            return i;
//        }
//    }
//    return types-1;
//}
//
//static double value(RGB rgb) {
//    return 255 - (rgb.R + rgb.G + rgb.B)/3;
//}
//
//void translate(const Bitmap& bmp, char const* out_name ) {
//    std::FILE* ascii = std::fopen(out_name, "w");
//
//    double val_err = 0.0;
//    for (int r=bmp.dims.height-1; r>=0; --r) {
//        for (int c=0; c<bmp.dims.width; ++c) {
//            val_err += value(bmp.data[r][c]);
//            int ind = index(val_err);
//            std::fputc(asciis[ind], ascii);
//            val_err -= thresholds[asciis[ind]-' '];
//        }
//        std::fputc('\n', ascii);
//    }
//
//    fclose(ascii);
//}
//
//void stretch(const Bitmap& bmp, Bitmap& bmp2) {
//    bmp2.allocate({bmp.dims.height/2, bmp.dims.width});
//    for (int r=0; r!=bmp.dims.height/2; ++r) {
//        for (int c=0; c!=bmp.dims.width; ++c) {
//            bmp2.data[r][c].R = (bmp.data[2*r][c].R + bmp.data[2*r+1][c].R)/2;
//            bmp2.data[r][c].G = (bmp.data[2*r][c].G + bmp.data[2*r+1][c].G)/2;
//            bmp2.data[r][c].B = (bmp.data[2*r][c].B + bmp.data[2*r+1][c].B)/2;
//        }
//    }
//}


