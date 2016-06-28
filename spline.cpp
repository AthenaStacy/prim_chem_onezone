#define integer int
#define  NH2DATA 41

void spline_coefficients(int nval, double *values, double (*coefficients)[NH2DATA-1]);
void spline_derivatives(int nval, double *values, double *derivatives);

/*=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////           S P L I N E _ E V A L           \\\\\\\\\\
c
=======================================================================
*/
 void   spline_eval(int nval, double *positions, double *values, int nnew,
                   double *new_positions, double *new_values)
{

    //void spline_coefficients(int nval, double *values, double *coefficients);
    
/*
c    written by: Simon Glover, January 2004
c
c  PURPOSE: Compute a spline fit to a set of nval data points.
c  Evaluate this fit at nnew specified locations and return the
c  nnew computed values.
*/
//      implicit NONE
//#include "cool.h"

    //integer nval, nnew;
    //REAL positions(nval), values(nval);
    //REAL new_positions(nnew), new_values(nnew);
    REAL coefficients[4][NH2DATA-1];
    REAL dpt, pt, a, b, c, d;
    integer I, J, index;

    
    
    spline_coefficients(nval, values, coefficients);
 
    for (I = 0; I<nnew; I++) {
        pt = new_positions[I];

//c We're using -1 to indicate that we've been asked for a value outside of
//c the range of our data, and setting the value for that point to zero;
//c anything fancier should probably be handled in the caller.
          
        if (pt < positions[0] || pt > positions[nval-1]) {
            index         = -1;
            new_values[I] = 0e0;
        }
        else if (pt == positions[nval-1]) {
//c Another special case, but this one is easy to deal with
            index         = -1;
            new_values[I] = values[nval-1];
        }
        else {
        index = 0;
            for(J = 0; J<nval-1; J++){
                if (pt >= positions[J] && pt < positions[J+1]) {
                    index = J;
                    dpt = (pt - positions[index]) / (positions[index+1] - positions[index]);
                    }
                }
        }

        if (index != -1) {
            a = coefficients[1-1][index];
            b = coefficients[2-1][index];
            c = coefficients[3-1][index];
            d = coefficients[4-1][index];
            new_values[I] = a + b * dpt + c * pow(dpt,2) + d * pow(dpt,3);
        }

     }

     return;
};
/*
c=======================================================================
c
c    \\\\\\\\\\        E  N  D   S U B R O U T I N E      //////////
c    //////////           S P L I N E _ E V A L           \\\\\\\\\\
c
c=======================================================================
*/
/*
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////   S P L I N E _ C O E F F I C I E N T S   \\\\\\\\\\
c
c=======================================================================
*/
    void spline_coefficients(int nval, double *values, double (*coefficients)[NH2DATA-1])
{
/*
c    written by: Simon Glover, January 2004
c
c  PURPOSE: Compute coefficients for a cubic spline fit to the set
c  of nval datapoints specified in values.
c
c  The solution follows the procedure of Bartels et al., 1998, "An 
c  Introduction to Splines for Use in Computer Graphics and Geometric 
c  Modelling", Morgan Kaufmann, ch. 3, pp 9-17, as outlined at 
c  http://mathworld.wolfram.com/CubicSpline.html
c
c  NB This routine does little or no error-checking, since it is only
c  intended to be used with a few small sets of known good data.
*/
//       implicit NONE
//#include "cool.h"
      //integer I, nval
      //REAL values(nval), derivatives(nval), coefficients(4,nval-1)
    //double coefficients[4][NH2DATA];
    int I;
    double derivatives[nval];
    
/*
c We first calculate the first derivative of the curve at each of
c our data points
*/
    
    spline_derivatives(nval, values, derivatives);

/*
c The spline coefficients can then be expressed in terms of the derivatives
c and the values at the data points
*/

    for( I = 0; I < nval-1; I++) {
        coefficients[1-1][I] = values[I];
        coefficients[2-1][I] = derivatives[I];
        coefficients[3-1][I] = 3.e0 * (values[I+1] - values[I]) - 2.e0 * derivatives[I] - derivatives[I+1];
        coefficients[4-1][I] = 2.e0 * (values[I] - values[I+1]) + derivatives[I] + derivatives[I+1];
       }

}
/*=======================================================================
c
c    \\\\\\\\\\        E N D     S U B R O U T I N E      //////////
c    //////////   S P L I N E _ C O E F F I C I E N T S   \\\\\\\\\\
c
c=======================================================================
*/

/*=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////    S P L I N E _ D E R I V A T I V E S    \\\\\\\\\\
c
c=======================================================================
*/
 void spline_derivatives(int nval, double *values, double *derivatives)
    {
/*
c    written by: Simon Glover, January 2004
c
c  PURPOSE: Compute the derivates of the spline curve required by 
c  the spline_coefficient subroutine. For more details, see the
c  preamble to that subroutine.
*/
//       implicit NONE
//#include "cool.h"

      //integer I, J, nval
      //REAL values(nval), derivatives(nval)
        int I, J;
        //double Cc[nval][nval+1];
        REAL f,sum;

        struct cc_data
        {
            double cc[NH2DATA+1];
        }*CC;
        
        CC=(struct cc_data *) malloc((nval)*sizeof(struct cc_data));
        
/* The derivatives we require are given by the solution of the matrix equation
c Ax = B, where:
c
c      |2 1        |       |D_1  |           |3*(y_1 - y_0)  |
c      |1 4 1      |       |D_2  |           |3*(y_2 - y_0)  |
c  A = |  .....    |,  x = |...  |  and  B = |3*(y_3 - y_1)  |
c      |    .....  |       |...  |           | ............. |
c      |      1 4 1|       |D_n-1|           |3*(y_n - y_n-2)|
c      |        1 2|       |D_n  |           |3*(y_n - y_n-1)|
c
c To solve this equation, we first construct the 'augmented', n by n+1 matrix
c Cc = AB:
*/
        for(I = 0; I < nval; I++){
            for(J = 0; J < nval+1; J++){
                CC[I].cc[J] = 0e0;
            }
        if (I == 0) {
            CC[I].cc[I]      = 2.e0;
            CC[I].cc[I+1]    = 1.e0;
            CC[I].cc[nval]   = 3.e0 * (values[2-1] - values[1-1]);
        }
        else if (I == nval) {
            CC[I].cc[I-1]    = 1.e0;
            CC[I].cc[I]      = 2.e0;
            CC[I].cc[nval]   = 3.e0 * (values[nval-1] - values[nval-1-1]);
        }
        else {
            CC[I].cc[I-1]    = 1.e0;
            CC[I].cc[I]      = 4.e0;
            CC[I].cc[I+1]    = 1.e0;
            CC[I].cc[nval]   = 3.e0 * (values[I+1] - values[I-1]);
        }
    }

// We then reduce this matrix to upper triangular form by Gaussian
//c elimination...
        for( I = 1; I < nval; I++) {
            f = CC[I].cc[I-1] / CC[I-1].cc[I-1];
            for( J = 0; J< nval+1; J++) {
                CC[I].cc[J] = CC[I].cc[J] - f * CC[I-1].cc[J];
            }
        }

//c ...and finally solve for the derivatives by back substitution
 //     derivatives(nval) = Cc(nval,nval+1) / Cc(nval,nval)
 
        for (I = nval-2; I >= 0; I--) {
            sum = 0e0;
            for( J = I; J < nval; J++){
                sum = sum + CC[I].cc[J] * derivatives[J];
                }
            derivatives[I] = (CC[I].cc[nval] - sum) / CC[I].cc[I];
            }

        free(CC);

    }
/*=======================================================================
c
c    \\\\\\\\\\        E N D     S U B R O U T I N E      //////////
c    //////////    S P L I N E _ D E R I V A T I V E S    \\\\\\\\\\
c
=======================================================================
*/
