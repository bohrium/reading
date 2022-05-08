/* author: samtenka
 * change: 2022-04-27 
 * create: 2018-03-29 
 * descrp: convert a given bitmap to an ASCII Art textfile
 * usage : compile (if needed) with `g++ -lm -std=c++11 Main.cpp -o main`,
           then run via `./main gbsmall.bmp` 
 */

#include "Bitmap.h"
#include "Ascii.h"
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) { std::cout << "Argument Error!"; return 0; }
    //init_thresholds("alpha.bmp", 8);

    Bitmap BMP, BMP2, BMP3;
    BMP.read_from(argv[1]);

    reframe         (BMP , BMP2);
    BMP2.write_to("stretched.bmp");

    clear_background(BMP2, BMP3, 0.15, 0.25);
    clear_background(BMP2, BMP3, 0.50, 1.5);
    //clear_background(BMP2, BMP3, 0.5, 2.5);

    BMP3.write_to("cleared.bmp");

    return 0;
}
