/* 
 c program:
 --------------------------------
  1. draws Mandelbrot set for Fc(z)=z*z +c
  using Mandelbrot algorithm ( boolean escape time )
 -------------------------------         
 2. technique of creating ppm file is  based on the code of Claudio Rocchini
 http://en.wikipedia.org/wiki/Image:Color_complex_plot.jpg
 create 24 bit color graphic file ,  portable pixmap file = PPM 
 see http://en.wikipedia.org/wiki/Portable_pixmap
 to see the file use external application ( graphic viewer)
  */
#include <stdio.h>
#include <math.h>
int main()
{
    /* screen ( integer) coordinate */
    int iX, iY;
    const int iXmax = 4000;
    const int iYmax = 4000;
    /* world ( double) coordinate = parameter plane*/
    double Cx, Cy;
    const double CxMin = -1.5;
    const double CxMax = 1.0;
    const double CyMin = -1.0;
    const double CyMax = 1.0;
    /* */
    double PixelWidth = (CxMax - CxMin) / iXmax;
    double PixelHeight = (CyMax - CyMin) / iYmax;
    /* color component ( R or G or B) is coded from 0 to 255 */
    /* it is 24 bit color RGB file */
    const int MaxColorComponentValue = 255;
    FILE *fp;
    char *filename = "new1.ppm";
    char *comment = "# "; /* comment should start with # */
    static unsigned char color[3];
    /* Z=Zx+Zy*i  ;   Z0 = 0 */
    double Zx, Zy;
    double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
    /*  */
    int Iteration;
    const int IterationMax = 500;
    /* bail-out value , radius of circle ;  */
    const double EscapeRadius = 100;
    double ER2 = EscapeRadius * EscapeRadius;
    /*create new file,give it a name and open it in binary mode  */
    fp = fopen(filename, "wb"); /* b -  binary mode */
    /*write ASCII header to the file*/
    fprintf(fp, "P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);
    /* compute and write image data bytes to the file*/
    for (iY = 0; iY < iYmax; iY++)
    {
        Cy = CyMin + iY * PixelHeight;
        if (fabs(Cy) < PixelHeight / 2)
            Cy = 0.0; /* Main antenna */
        for (iX = 0; iX < iXmax; iX++)
        {
            Cx = CxMin + iX * PixelWidth;
            /* initial value of orbit = critical point Z= 0 */
            Zx = 0.0;
            Zy = 0.0;
            Zx2 = Zx * Zx;
            Zy2 = Zy * Zy;
            /* */
            for (Iteration = 0; Iteration < IterationMax && ((Zx2 + Zy2) < ER2); Iteration++)
            {
                Zy = 2 * Zx * Zy + Cy;
                Zx = Zx2 - Zy2 + Cx;
                Zx2 = Zx * Zx;
                Zy2 = Zy * Zy;
            };

            if (Iteration == IterationMax)
            {
                // Point within the set. Mark it as black
                color[0] = 0;
                color[1] = 0;
                color[2] = 0;
            }
            else
            {
                double c = 3 * log((double)Iteration) / log((double)(IterationMax)-1.0);
                if (c < 1)
                {
                    color[0] = 0;
                    color[1] = 0;
                    color[2] = 255 * c;
                }
                else if (c < 2)
                {
                    color[0] = 0;
                    color[1] = 255 * (c - 1);
                    color[2] = 255;
                }
                else
                {
                    color[0] = 255 * (c - 2);
                    color[1] = 255;
                    color[2] = 255;
                }
            }
            /*write color to the file*/
            fwrite(color, 1, 3, fp);
        }
    }
    fclose(fp);
    return 0;
}