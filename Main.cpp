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
    load_fonts("fonts.bmp");

    Bitmap BMP, BMP2, BMP3;
    BMP.read_from(argv[1]);

    ocr(BMP);

    //reframe         (BMP , BMP2);
    //BMP2.write_to("stretched.bmp");

    //clear_background(BMP2, BMP3, 0.5 , 1.5 );
    //BMP3.write_to("cleared.bmp");

    //ocr(BMP3);

    return 0;
}
