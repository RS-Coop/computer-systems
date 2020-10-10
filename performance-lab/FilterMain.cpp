#include <stdio.h>
#include "cs1300bmp.h"
#include <iostream>
#include <fstream>
#include "Filter.h"

// #include <opencv2/opencv.hpp> //For openCV

using namespace std;

#include "rdtsc.h"

//
// Forward declare the functions
//
Filter * readFilter(string filename);
double applyFilter(Filter *filter, cs1300bmp *input, cs1300bmp *output);

int
main(int argc, char **argv)
{

  if ( argc < 2) {
    fprintf(stderr,"Usage: %s filter inputfile1 inputfile2 .... \n", argv[0]);
  }

  //
  // Convert to C++ strings to simplify manipulation
  //
  string filtername = argv[1];

  //
  // remove any ".filter" in the filtername
  //
  string filterOutputName = filtername;
  string::size_type loc = filterOutputName.find(".filter");
  if (loc != string::npos) {
    //
    // Remove the ".filter" name, which should occur on all the provided filters
    //
    filterOutputName = filtername.substr(0, loc);
  }

  Filter *filter = readFilter(filtername);

  double sum = 0.0;
  int samples = 0;

  for (int inNum = 2; inNum < argc; inNum++) {
    string inputFilename = argv[inNum];
    string outputFilename = "filtered-" + filterOutputName + "-" + inputFilename;
    struct cs1300bmp *input = new struct cs1300bmp;
    struct cs1300bmp *output = new struct cs1300bmp;
    int ok = cs1300bmp_readfile( (char *) inputFilename.c_str(), input);

    if ( ok ) {
      double sample = applyFilter(filter, input, output);
      sum += sample;
      samples++;
      cs1300bmp_writefile((char *) outputFilename.c_str(), output);
    }
    delete input;
    delete output;
  }
  fprintf(stdout, "Average cycles per sample is %f\n", sum / samples);

}

class Filter *
readFilter(string filename)
{
  ifstream input(filename.c_str());

  if ( ! input.bad() ) {
    int size = 0;
    input >> size;
    Filter *filter = new Filter(size);
    int div;
    input >> div;
    filter -> setDivisor(div);
    for (int i=0; i < size; i++) {
      for (int j=0; j < size; j++) {
	int value;
	input >> value;
	filter -> set(i,j,value);
      }
    }
    return filter;
  } else {
    cerr << "Bad input in readFilter:" << filename << endl;
    exit(-1);
  }
}


double
applyFilter(class Filter *filter, cs1300bmp *input, cs1300bmp *output)
{

  long long cycStart, cycStop;

  cycStart = rdtscll();
    
  //Some other speed up vars
    int ih = input -> height; //Replaces access in first two for loops
    int iw = input -> width; //Replaces access in first two for loops

  output -> width = iw;
  output -> height = ih;
    
  //Replaces upper bound in for loop  
    int iw_2 = iw - 2;
    int ih_2 = ih - 2;
    
  //Local variables for speed up
//     int filter_size = filter -> getSize(); //Replaces call in innermost two forloops
//     int filter_divisor = filter -> getDivisor(); //Replaces call before first if
//     int* filter_data = filter -> getArr(); //Replaces innermost for loop call
//     int dims = filter -> getSize();
//     int filter_size = filter -> dim;
    int filter_divisor = filter -> divisor;
    int* filter_data = filter -> data;
    
  //Filter data increment
    short n; //
    int dims_square = 8192 * 8192;
    
  //Pointers to output and input color arrays
//     int(*output_arr)[8192][8192] = output -> color; //Replaces all output -> color[plane][row][col] references
//     int(*input_arr)[8192][8192] = input -> color;
//     int* input_planes;
    int* output_flat = output -> color_flat;
    int* input_flat_p1 = input -> color_flat;
    
  //Changing these two specific plane pointers.
    int* input_flat_p2 = input_flat_p1 + dims_square;
    int* input_flat_p3 = input_flat_p2 + dims_square;
    
  //Local register stored sum maybe.
    int conv2D_p1; //Replaces most output_arr[][][] references.
    int conv2D_p2;
    int conv2D_p3;
    
  //Reduce add ops
    int ri;
    
  //Reduce memory access for filter data
    short f0;
    short f1;
    short f2;
    
  #pragma omp parallel for
  for(int row = 0; row < ih_2; row++) {
    for(int col = 0; col < iw_2 ; col++) {

	conv2D_p1 = 0;
    conv2D_p2 = 0;
    conv2D_p3 = 0;
    n = 0;
    ri = 8192*row ;
        
    for (short i = 0; i < 3; i++) { //COLOR_RED + 3*row + 8192*8192*col
//         ri = 8192*(row + i);
        f0 = *(filter_data + n);
        f1 = *(filter_data + n + 1);
        f2 = *(filter_data + n + 2);
        
        //Flat array stuff
        
//         conv2D_p1 += (*(input_flat_p1 + ri + col) * f0)
//             + (*(input_flat_p1 + ri + col+1) * f1)
//             + (*(input_flat_p1 + ri + col+2) * f2);
        
//         conv2D_p2 += (*(input_flat_p2 + ri + col) * f0)
//             + (*(input_flat_p2 + ri + col+1) * f1)
//             + (*(input_flat_p2 + ri + col+2) * f2);
        
//         conv2D_p3 += (*(input_flat_p3 + ri + col) * f0)
//             + (*(input_flat_p3 + ri + col+1) * f1)
//             + (*(input_flat_p3 + ri + col+2) * f2);
        
        conv2D_p1 += (*(input_flat_p1 + ri + col) * f0)
            + (*(input_flat_p1 + ri + col+1) * f1)
            + (*(input_flat_p1 + ri + col+2) * f2);
        
        conv2D_p2 += (*(input_flat_p2 + ri + col) * f0)
            + (*(input_flat_p2 + ri + col+1) * f1)
            + (*(input_flat_p2 + ri + col+2) * f2);
        
        conv2D_p3 += (*(input_flat_p3 + ri + col) * f0)
            + (*(input_flat_p3 + ri + col+1) * f1)
            + (*(input_flat_p3 + ri + col+2) * f2);

        n += 3;
        
        ri += 8192;
	}
	
    if(filter_divisor != 1)
    {
        conv2D_p1 /= filter_divisor;
        conv2D_p2 /= filter_divisor;
        conv2D_p3 /= filter_divisor;
    }
          
    conv2D_p1 = conv2D_p1 < 0 ? 0 : conv2D_p1;
    conv2D_p1 = conv2D_p1 > 255 ? 255 : conv2D_p1;
          
    conv2D_p2 = conv2D_p2 < 0 ? 0 : conv2D_p2;
    conv2D_p2 = conv2D_p2 > 255 ? 255 : conv2D_p2;
          
    conv2D_p3 = conv2D_p3 < 0 ? 0 : conv2D_p3;
    conv2D_p3 = conv2D_p3 > 255 ? 255 : conv2D_p3;

        
        //Flat array stuff
        *(output_flat + 8192*(row+1) + col+1) = conv2D_p1;
        *(output_flat + dims_square + 8192*(row+1) + col+1) = conv2D_p2;
        *(output_flat + 2*dims_square + 8192*(row+1) + col+1) = conv2D_p3;

    }
      
//       input_flat_p1 += 8192;
//       input_flat_p2 += 8192;
//       input_flat_p3 += 8192;
  }

  cycStop = rdtscll();
  double diff = cycStop - cycStart;
  double diffPerPixel = diff / (output -> width * output -> height);
  fprintf(stderr, "Took %f cycles to process, or %f cycles per pixel\n",
	  diff, diff / (output -> width * output -> height));
  return diffPerPixel;
}
