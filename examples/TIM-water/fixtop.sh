#!/bin/bash

file=$1
sed -i 's/      8    glu_C8      1     GLU    CG      3  -0.220000     12.011   glh_C8  -0.120000     12.011/     8   opls_274      1    GLU     CG      3  -0.220000     12.011  opls_136  -0.120000     12.011/' $file

sed -i 's/      9    glu_H9      1     GLU   HG1      3   0.060000      1.008   glh_H9   0.060000      1.008/     9   opls_140      1    GLU    HG1      3   0.060000      1.008  opls_140   0.060000      1.008/' $file

sed -i 's/     10   glu_H10      1     GLU   HG2      3   0.060000      1.008  glh_H10   0.060000      1.008/    10   opls_140      1    GLU    HG2      3   0.060000      1.008  opls_140   0.060000      1.008/' $file

sed -i 's/     11   glu_C11      1     GLU    CD      4   0.700000     12.011  glh_C11   0.520000     12.011/    11   opls_271      1    GLU     CD      4   0.700000     12.011  opls_267   0.520000     12.011/' $file

sed -i 's/     12   glu_O12      1     GLU   OE1      4  -0.800000    15.9994  glh_O12  -0.440000    15.9994/    12   opls_272      1    GLU    OE1      4  -0.800000    15.9994  opls_269  -0.440000    15.9994/' $file

sed -i 's/     13   glu_O13      1     GLU   OE2      4  -0.800000    15.9994  glh_O13  -0.530000    15.9994/    13   opls_272      1    GLU    OE2      4  -0.800000    15.9994  opls_268  -0.530000    15.9994/' $file

sed -i 's/     30   sub_H15      2     SUB   H15      6   0.060000      1.008  glh_H14   0.450000      1.008   ; qtot -3/     30   sub_H15      2     SUB   H15      6   0.060000      1.008 opls_270   0.450000      1.008   ; qtot -3/' $file


#sed -i 's///' $file


