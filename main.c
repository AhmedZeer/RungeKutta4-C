#include <stdio.h>
#include <stdlib.h>

////////////////////// SOLVING ODES USING RUNGE KUTTA-4 ///////////////////// 


/* This Struct will act as our polynom 
 * coef stands for coefficient.
 * exp stands for the biggest exponent. */
typedef struct {
    double *coef;
    int exp;
} Polynomial;


/*Evaluate our polynom at a specific point using Horner's Rule*/
double EvalPoli( const Polynomial *poli , double x)  
{
    double result = 0.0;
    int i;

    for ( i = poli->exp ; i >= 0; i--)              
    {
        result = result*x + poli->coef[i];
    }
    
    return result;
}

/*Parameters Struct*/
typedef struct{
    double x0;                                                          
    double y0;                                                         
    double stepSize;                                                    
    int numOfSteps;  
    double (*dydx)( double, const Polynomial *poli, double, double);   
} Param; 

/*This will act as our diffrential equation dy/dx*/
double dydx( double a, const Polynomial *poli, double x, double y)
{
    return (a*y + EvalPoli(poli , x));
}

double rungeKutta( double a, const Polynomial *poli, Param param ) { 
    double x = param.x0;
    double y = param.y0;
    int i;

    for ( i = 0; i < param.numOfSteps; ++i) {
        double k1 = param.stepSize * param.dydx( a , poli ,x, y);
        double k2 = param.stepSize * param.dydx( a , poli ,x + 0.5 * param.stepSize, y + 0.5 * k1);
        double k3 = param.stepSize * param.dydx( a , poli ,x + 0.5 * param.stepSize, y + 0.5 * k2);
        double k4 = param.stepSize * param.dydx( a , poli ,x + param.stepSize, y + k3);

        printf("\n\ny: %lf", y);
        printf("\nx: %lf", x);

        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        x += param.stepSize;
    }

    return y;
}

int main()
{
    printf("\n NOTE: your diff equation shall be like this y' = A*y + p(x).\n");
    printf("\nP(X) will act as out Polynomial.\n");

    Polynomial poli;
    double x0;
    double y0;
    double x;
    double a;
    int exp, i;
    double stepSize = 0.1;
    double numOfSteps;

    printf("\nPolynomial's Degree: ");
    scanf("%d", &exp);

    poli.coef = (double*)malloc((exp)*sizeof(int));
    if( poli.coef == NULL)
    {
        printf("MEMORY ALLOCATION FAILED !");
        return 1;
    }

    poli.exp = exp;

    printf("\nConstant: ");
    scanf("%lf", &poli.coef[0]);

    for ( i = 1; i <= exp; i++)
    {
        printf("\na * x^%d term's coefficient: ", i);
        scanf("%lf", &poli.coef[i]);
    }
    
    printf("\nA * y + p(x) equation's A value: ");
    scanf("%lf", &a);

    printf("\nx0: ");
    scanf("%lf", &x0);

    printf("\ny0: ");
    scanf("%lf", &y0);
    
    printf("\nx value you want to find: ");
    scanf("%lf", &x);

    numOfSteps = ((x - x0) / stepSize);

    Param param = {
        .dydx = dydx,
        .x0 = x0,
        .y0 = y0,
        .stepSize = stepSize,
        .numOfSteps = numOfSteps,
    };

    double result = rungeKutta(a, &poli , param);
    printf("\n====================================");
    printf("\nResult: %lf\n", result);
    
    system("pause");

    return 0;
}
