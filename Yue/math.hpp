//
//  math.hpp
//  CGNew2
//
//  Created by 吴越 on 16/3/11.
//  Copyright © 2016年 Yue Wu. All rights reserved.
//

#ifndef math_hpp
#define math_hpp

#include <stdio.h>
class SparseMatrix;

class Math {
public:
    static void conjugate_gradient(
                                   const SparseMatrix &a, const double *b, int num_iterations, double *x);
};
#endif /* math_hpp */
