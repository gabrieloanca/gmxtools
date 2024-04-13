#!/bin/bash

file=$1
sed -i 's/   2544    glu_C8    164     GLU    CG    858  -0.220000     12.011   glh_C8  -0.120000     12.011/  2544   opls_274    164    GLU     CG    858  -0.220000     12.011  opls_136  -0.120000     12.011/' $file

sed -i 's/   2545    glu_H9    164     GLU   HG1    858   0.060000      1.008   glh_H9   0.060000      1.008/  2545   opls_140    164    GLU    HG1    858   0.060000      1.008  opls_140   0.060000      1.008/' $file

sed -i 's/   2546   glu_H10    164     GLU   HG2    858   0.060000      1.008  glh_H10   0.060000      1.008/  2546   opls_140    164    GLU    HG2    858   0.060000      1.008  opls_140   0.060000      1.008/' $file

sed -i 's/   2547   glu_C11    164     GLU    CD    859   0.700000     12.011  glh_C11   0.520000     12.011/  2547   opls_271    164    GLU     CD    859   0.700000     12.011  opls_267   0.520000     12.011/' $file

sed -i 's/   2548   glu_O12    164     GLU   OE1    859  -0.800000    15.9994  glh_O12  -0.440000    15.9994/  2548   opls_272    164    GLU    OE1    859  -0.800000    15.9994  opls_269  -0.440000    15.9994/' $file

sed -i 's/   2549   glu_O13    164     GLU   OE2    859  -0.800000    15.9994  glh_O13  -0.530000    15.9994/  2549   opls_272    164    GLU    OE2    859  -0.800000    15.9994  opls_268  -0.530000    15.9994/' $file

sed -i 's/   7557   sub_H15    495     SUB   H15   2543   0.060000      1.008  glh_H14   0.450000      1.008   ; qtot -8/   7557   sub_H15    495     SUB   H15   2543   0.060000      1.008 opls_270   0.450000      1.008   ; qtot -8/' $file

#sed -i 's///' $file

