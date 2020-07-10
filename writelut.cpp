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

// [[Rcpp::export]]
int writelut(IntegerVector x, std::string filename) {
   // Rcout << "Lorem" << std::endl;
    // std::ofstream ofs ("/Users/puetz/test.txt", std::ofstream::out);
    // 
    // ofs << "lorem ipsum" << std::endl;
    // 
    // ofs.close();
    
    //Rcout << filename << std::endl;
    std::ofstream out(filename, 
                      std::ofstream::out | std::ofstream::binary);
    int count = 0;
    for(int i:x){
        if(i > 255){
            Rcerr << "element > 255 encountered: "<< i << std::endl;
            out.close();
            return 99;
        }
        //if(count++ < 6){
            Rcout << i << ' ';
            if(++count%3 == 0)
                Rcout << "  ";
        //}
        out << (unsigned char)i;
            
    }
    out.close();
    Rcout << std::endl;
    return 0;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
writelut(as.integer(sample(0:255, 18)), 'min.lut')
*/
