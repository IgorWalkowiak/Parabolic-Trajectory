#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <math.h>

const double G=9.80665;
double B=1.0;

int funcX (double t, const double y[], double f[], void *params)
{
  (void)(t);
  double mu=*((double *)params+0);

  double Ux=*((double *)params+1);
  double Uy=*((double *)params+2);
  double Wx=*((double *)params+3);
  double Wy=*((double *)params+4);

  f[0] = y[1];
  f[1] = -mu*sqrt((y[1]-Wx)*(y[1]-Wx)+(Uy-Wy)*(Uy-Wy))*(y[1]-Wx);

  return GSL_SUCCESS;
}

int jacX (double t, const double y[], double *dfdy,
     double dfdt[], void *params)
{
  (void)(t);
    double mu=*((double *)params+0);
    double Ux=*((double *)params+1);
    double Uy=*((double *)params+2);
    double Wx=*((double *)params+3);
    double Wy=*((double *)params+4);
    double x=*((double *)params+5);
    //double x=-mu*(Ux-Wx)*(Ux-Wx)/sqrt((Ux-Wx)(Ux-Wx)+(Uy-Wy)*(Uy-Wy))-mu*sqrt((Ux-Wx)*(Ux-Wx)+(Uy-Wy)*(Uy-Wy))
  gsl_matrix_view dfdy_mat
    = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, -mu*(Ux-Wx)*(Ux-Wx)/sqrt((Ux-Wx)*(Ux-Wx)+(Uy-Wy)*(Uy-Wy))-mu*sqrt((Ux-Wx)*(Ux-Wx)+(Uy-Wy)*(Uy-Wy)));
  gsl_matrix_set (m, 1, 1, 0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

//___________Y_______________

int funcY (double t, const double y[], double f[], void *params)
{
  (void)(t);

    double mu=*((double *)params+0);
    double Ux=*((double *)params+1);
    double Uy=*((double *)params+2);
    double Wx=*((double *)params+3);
    double Wy=*((double *)params+4);

  f[0] = y[1];
  f[1] = -G-mu*sqrt((Ux-Wx)*(Ux-Wx)+(y[1]-Wy)*(y[1]-Wy))*(y[1]-Wy);
  return GSL_SUCCESS;
}

int jacY (double t, const double y[], double *dfdy,
     double dfdt[], void *params)
{
  (void)(t);
    double mu = *(double *)params;
    double Ux=*((double *)params+1);
    double Uy=*((double *)params+2);
    double Wx=*((double *)params+3);
    double Wy=*((double *)params+4);
    double x=*((double *)params+5);
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, -mu*(Uy-Wy)*(Uy-Wy)/sqrt((Ux-Wx)*(Ux-Wx)+(Uy-Wy)*(Uy-Wy))-mu*sqrt((Ux-Wx)*(Ux-Wx)+(Uy-Wy)*(Uy-Wy)));
  gsl_matrix_set (m, 1, 1, 0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}



int main (void)
{

    double Ukat=25;
    double Upredkosc=30;
    double Wkat=170;
    double Wpredkosc=3;
   /*   std::cout<<"Podaj predkosc poczatkowa i kat w stopniach";
        std::cin>>Upredkosc>>Ukat;

        std::cout<<"Podaj predkosc wiatru i kat w stopniach";
        std::cin>>Wpredkosc>>Wkat;
    */


    Ukat=Ukat*0.017453;
    Wkat=Wkat*0.017453;

    double k=0.1;

    double Ux=Upredkosc*cos(Ukat);
    double Uy=Upredkosc*sin(Ukat);
    double Wx=Wpredkosc*cos(Wkat);
    double Wy=Wpredkosc*sin(Wkat);

    double Fx=-k*sqrt((Ux-Wx)*(Ux-Wx)+(Uy-Wy)*(Uy-Wy))*(Ux-Wx);
    double Fy=-k*sqrt((Ux-Wx)*(Ux-Wx)+(Uy-Wy)*(Uy-Wy))*(Uy-Wy);



    double par[6] = {k, Ux, Uy, Wx, Wy,1};
    //par[5]=-k*(Ux-Wx)*(Ux-Wx)/sqrt((Ux-Wx)*(Ux-Wx)+(Uy-Wy)*(Uy-Wy))-k*sqrt((Ux-Wx)*(Ux-Wx)+(Uy-Wy)*(Uy-Wy));
    //par[6]=-k*(Uy-Wy)*(Uy-Wy)/sqrt((Ux-Wx)*(Ux-Wx)+(Uy-Wy)*(Uy-Wy))-k*sqrt((Ux-Wx)*(Ux-Wx)+(Uy-Wy)*(Uy-Wy));
    double t = 0.0, t1 = 10.0;
    double x[2] = { 0.0, Upredkosc*cos(Ukat) };
    double y[2] = { 0.0, Upredkosc*sin(Ukat) };
    int j;

        gsl_odeiv2_system sysX = {funcX, jacX, 2, &par};
        gsl_odeiv2_driver * X = gsl_odeiv2_driver_alloc_y_new(&sysX, gsl_odeiv2_step_rk8pd,1e-6, 1e-6, 0.0);

        gsl_odeiv2_system sysY = {funcY, jacY, 2, &par};
        gsl_odeiv2_driver * Y = gsl_odeiv2_driver_alloc_y_new(&sysY, gsl_odeiv2_step_rk8pd,1e-6, 1e-6, 0.0);

  for (int i = 1; i <= 10000; i++)
    {

        int s= gsl_odeiv2_driver_apply_fixed_step(X,&t,1e-5,1000,x);
        par[1]=x[1];
        par[5]=-par[0]*(par[1]-par[3])*(par[1]-par[3])/sqrt((par[1]-par[3])*(par[1]-par[3])+(par[2]-par[4])*(par[2]-par[4]))-par[0]*sqrt((par[1]-par[3])*(par[1]-par[3])+(par[2]-par[4])*(par[2]-par[4]));
        par[6]=-par[0]*(par[2]-par[4])*(par[2]-par[4])/sqrt((par[1]-par[3])*(par[1]-par[3])+(par[2]-par[4])*(par[2]-par[4]))-par[0]*sqrt((par[1]-par[3])*(par[1]-par[3])+(par[2]-par[4])*(par[2]-par[4]));
        if (s != GSL_SUCCESS)
        {
          printf ("error: driver returned %d\n", s);
          break;
        }

        s= gsl_odeiv2_driver_apply_fixed_step(Y,&t,1e-5,1000,y);
        par[2]=y[1];
        par[5]=-par[0]*(par[1]-par[3])*(par[1]-par[3])/sqrt((par[1]-par[3])*(par[1]-par[3])+(par[2]-par[4])*(par[2]-par[4]))-par[0]*sqrt((par[1]-par[3])*(par[1]-par[3])+(par[2]-par[4])*(par[2]-par[4]));
        par[6]=-par[0]*(par[2]-par[4])*(par[2]-par[4])/sqrt((par[1]-par[3])*(par[1]-par[3])+(par[2]-par[4])*(par[2]-par[4]))-par[0]*sqrt((par[1]-par[3])*(par[1]-par[3])+(par[2]-par[4])*(par[2]-par[4]));
        if (s != GSL_SUCCESS)
        {
            printf ("error: driver returned %d\n", s);
            break;
        }


    std::cout<<t<<" "<<y[0]<<" "<<y[1]<<" "<<x[0]<<" "<<x[1]<<"\n";
    if(y[0]<0){j=i;break;}
    }
   gsl_odeiv2_driver_free (Y);

  return 0;
}



/*
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <math.h>

const double G=9.80665;
double B=1.0;

int funcX (double t, const double y[], double f[], void *params)
{
  (void)(t);
  double mu = *(double *)params;
  f[0] = y[1];
  f[1] = -mu*y[1];
  return GSL_SUCCESS;
}

int jacX (double t, const double y[], double *dfdy,
     double dfdt[], void *params)
{
  (void)(t);
  double mu = *(double *)params;
  gsl_matrix_view dfdy_mat
    = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, -mu);
  gsl_matrix_set (m, 1, 1, 0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

//___________Y_______________

int funcY (double t, const double y[], double f[], void *params)
{
  (void)(t);
  double mu = *(double *)params;
  f[0] = y[1];
  f[1] = -G-mu*y[1];
  return GSL_SUCCESS;
}

int jacY (double t, const double y[], double *dfdy,
     double dfdt[], void *params)
{
  (void)(t);
  double mu = *(double *)params;
  gsl_matrix_view dfdy_mat
    = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, -mu);
  gsl_matrix_set (m, 1, 1, 0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}



int main (void)
{

    double kat=1.1;
    double predkosc=50;

    double mu = 0.5;

    double t = 0.0, t1 = 10.0;
    double x[2] = { 0.0, predkosc*cos(kat) };
    double y[2] = { 0.0, predkosc*sin(kat) };
    int j;

        gsl_odeiv2_system sysX = {funcX, jacX, 2, &mu};
        gsl_odeiv2_driver * X = gsl_odeiv2_driver_alloc_y_new(&sysX, gsl_odeiv2_step_rk8pd,1e-6, 1e-6, 0.0);

        gsl_odeiv2_system sysY = {funcY, jacY, 2, &mu};
        gsl_odeiv2_driver * Y = gsl_odeiv2_driver_alloc_y_new(&sysY, gsl_odeiv2_step_rk8pd,1e-6, 1e-6, 0.0);

  for (int i = 1; i <= 1000; i++)
    {

        int s= gsl_odeiv2_driver_apply_fixed_step(X,&t,1e-5,1000,x);
          s= gsl_odeiv2_driver_apply_fixed_step(Y,&t,1e-5,1000,y);


         gsl_odeiv2_driver_reset(Y);




    std::cout<<t<<" "<<y[0]<<" "<<y[1]<<" "<<x[0]<<" "<<x[1]<<"\n";
    if(y[0]<0){j=i;break;}
    }
   gsl_odeiv2_driver_free (Y);

  return 0;
}


*/




