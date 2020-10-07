#include <Rcpp.h>
using namespace Rcpp;

#include<iostream>
#include<iomanip>
#include<fstream>

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' write LUT
//' 
//' Write RGB color lookup table (LUT) of 3\eqn{\times}{x}256 bytes. 
//' Implemented with Rcpp to allow writing bytes.
//' 
//' @name writelut
//' @param x vector of 768 bytes (treated as unsigned integer)
//' @param filename name for output file
//' 
//' @return 0 on success
//' @export
// [[Rcpp::export]]
int writelut(IntegerVector x, std::string filename) {
   
    //Rcout << filename << std::endl;
    std::ofstream out(filename, 
                      std::ofstream::out | std::ofstream::binary);
    if(out.fail())    //bad() function will check for badbit
      {
        //Rcerr << "Opening file " << filename << " failed" << std::endl;
        return 1;
      }
    // int count = 0;
    for(int i:x){
        if(i > 255){
            //Rcerr << "element > 255 encountered: "<< i << std::endl;
            out.close();
            // remove corrupt LUT file
            std::string cmd("rm ");
            cmd.append(filename);
            const char *command = cmd.c_str(); 
            std::system(command);
            return 99;
        }   
        //if(count++ < 6){
            // Rcout << i << ' ';
            // if(++count%3 == 0)
            //     Rcout << "  ";
        //}
        out << (unsigned char)i;
            
    }
    out.close();
    if(out.bad())    //bad() function will check for badbit
    {
        //Rcerr << "Error writing to " << filename << std::endl;
        return 2;
    }
    // Rcout << std::endl;
    return 0;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
writelut(as.integer(sample(0:255, 18)), 'min.lut')
*/

