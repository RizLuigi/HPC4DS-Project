#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
//#include <omp.h>
#include <math.h>

int main(int argc, char *argv[])
{
    int NUM_ZOOMS = strtod(argv[1], NULL);
    int NUM_THREADS = strtod(argv[2], NULL);

    int my_rank, comm_sz;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (comm_sz <= NUM_ZOOMS)
    {
        // In questo caso i processi sono indipendenti: ognuno lavora su uno o più zoom
        for (int i = 0; i < NUM_ZOOMS / comm_sz; i++)
        {
            int actual_zoom = my_rank * NUM_ZOOMS / comm_sz + i;
            printf("Process %d working on zoom %d\n", my_rank, actual_zoom);
        }
    }
    else
    {
        // In questo caso più processi lavorano su uno stesso zoom, quindi ognuno avrà un master di riferimento
        // Forse è il caso  di creare nuovi communicators in modo da poter utilizzare collective communications?
        // https://mpitutorial.com/tutorials/introduction-to-groups-and-communicators/
        int proc_per_zoom = comm_sz / NUM_ZOOMS;
        int actual_zoom = my_rank / proc_per_zoom;
        int my_master = my_rank - my_rank % proc_per_zoom;
        printf("Process %d working on zoom %d with master %d\n", my_rank, actual_zoom, my_master);

        /* screen ( integer) coordinate */
        int iX, iY;
        const int iXmax = 40;
        const int iYmax = 40;
        /* world ( double) coordinate = parameter plane*/
        double Cx, Cy;
        const double CxMin = -1.0;
        const double CxMax = 1.0;
        const double CyMin = -1.0;
        const double CyMax = 1.0;
        /* */
        double PixelWidth = (CxMax - CxMin) / iXmax;
        double PixelHeight = (CyMax - CyMin) / iYmax;
        /* color component ( R or G or B) is coded from 0 to 255 */
        /* it is 24 bit color RGB file */
        const int MaxColorComponentValue = 255;
        
        //FILE *fp;
        //char *filename = "new1.ppm";
        //char *comment = "# "; /* comment should start with # */
        static unsigned char color[3];  //OCCHIO ALLO STATIC
        unsigned char ppmMatrix[iXmax][iYmax][3];


        /* Z=Zx+Zy*i  ;   Z0 = 0 */
        double Zx, Zy;
        double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
        /*  */
        int Iteration;
        const int IterationMax = 50;
        /* bail-out value , radius of circle ;  */
        const double EscapeRadius = 100;
        double ER2 = EscapeRadius * EscapeRadius;
        /*create new file,give it a name and open it in binary mode  */
        //fp = fopen(filename, "wb"); /* b -  binary mode */
        /*write ASCII header to the file*/
        //fprintf(fp, "P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);
        /* compute and write image data bytes to the file*/
        printf("INIZIO FOR");
        for (iY = my_rank; iY < iYmax; iY+comm_sz)
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
                //fwrite(color, 1, 3, fp);
                for (int i=0; i<3; i++) {
                    ppmMatrix[iY][iX][i] = color[i];
                }
            }
        }
        printf("FINE FOR\n");
        //fclose(fp);
        for(int iY=my_rank; iY<iYmax; iY+comm_sz) {
            for(int iX=0; iX<iXmax; iX++) {
                for(int i=0; i<3; i++) {
                    printf("%c,", ppmMatrix[iX][iY][i]);
                }
                printf("\t");
            }
            printf("\n");
        }
    }

    MPI_Finalize();
    return 0;
}
