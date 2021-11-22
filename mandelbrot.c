#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
//#include <omp.h>
#include <math.h>
#include <argp.h>

struct argp_option options[] =
{
    {"threads", 't', "threadCount", 0, "Use threadCount threads per process"},
    {"zooms", 'z', "zoomCount", 0 , "Create zoomCount images"},
    {"zoom-factor", 'f', "zoomFactor", 0, "Use zoomFactor as magnification factor"},
    {"final-x", 'x', "X", 0, "Use X as final x-coordinate"},
    {"final-y", 'y', "Y", 0, "Use Y as final y-coordinate"},
    {0 }
};

struct arguments
{
  int threads, zoom, zoom_factor;
  double final_x, final_y;
};

static int parse_opt (int key, char* arg, struct argp_state *state)
{   
    struct arguments *arguments = state->input;

    switch (key)
    {
    case 't':
        arguments->threads = atoi(arg);
        break;
    case 'z':
        arguments->zoom = atoi(arg);
        break;
    case 'f':
        arguments->zoom_factor = atoi(arg);
        break;
    case 'x':
        arguments->final_x = atoi(arg);
        break;
    case 'y':
        arguments->final_y = atoi(arg);
        break;
    default:
      return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt};

int main(int argc, char *argv[])
{
    struct arguments arguments;

    // DEFAULTS
    arguments.threads = 10;
    arguments.zoom = 4;
    arguments.zoom_factor = 5;
    arguments.final_x = -0.7;
    arguments.final_y = 0.26;

    argp_parse (&argp, argc, argv, 0, 0, &arguments);

    int my_rank, comm_sz;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // parametri per lo zoom
    double final_x = arguments.final_x;
    double final_y = arguments.final_y;
    double zoom_inc = arguments.zoom_factor;

    // default values
    /* screen ( integer) coordinate */
    int iX, iY;
    const int iXmax = 800;
    const int iYmax = 800;
    /* world ( double) coordinate = parameter plane*/
    double Cx, Cy;
    const double CxMin = -2.5;
    const double CxMax = 1.5;
    const double CyMin = -2.0;
    const double CyMax = 2.0;

    if (comm_sz <= arguments.zoom)
    {
        // In questo caso i processi sono indipendenti: ognuno lavora su uno o più zoom
        for (int i = 0; i < arguments.zoom / comm_sz; i++)
        {
            int actual_zoom = my_rank * arguments.zoom / comm_sz + i;
            printf("Process %d working on zoom %d\n", my_rank, actual_zoom);
        }
    }
    else
    {
        // In questo caso più processi lavorano su uno stesso zoom, quindi ognuno avrà un master di riferimento
        // Forse è il caso  di creare nuovi communicators in modo da poter utilizzare collective communications?
        // https://mpitutorial.com/tutorials/introduction-to-groups-and-communicators/
        int proc_per_zoom = comm_sz / arguments.zoom;
        int actual_zoom = my_rank / proc_per_zoom;
        int my_master = my_rank - my_rank % proc_per_zoom;
        printf("Process %d working on zoom %d with master %d\n", my_rank, actual_zoom, my_master);

        // creazione communicators gropus
        int working_frame = my_rank / proc_per_zoom;
        MPI_Comm frame_comm;
        MPI_Comm_split(MPI_COMM_WORLD, working_frame, comm_sz, &frame_comm);
        int frame_rank, frame_size;
        MPI_Comm_rank(frame_comm, &frame_rank);
        MPI_Comm_size(frame_comm, &frame_size);
        int frame_maser = frame_rank - frame_rank % proc_per_zoom;

        printf("%d/%d - frame_size: %d - frame_rank: %d frame_maser: %d\n", my_rank, actual_zoom, frame_size, frame_rank, frame_maser);

        // new coordinates for current frame
        double CxMax_cur = CxMax;
        double CxMin_cur = CxMin;
        double CyMax_cur = CyMax;
        double CyMin_cur = CyMin;

        if (actual_zoom != 0)
        {
            int cur_zoom = zoom_inc * actual_zoom;
            CxMax_cur = CxMax_cur / cur_zoom + final_x;
            CyMax_cur = CyMax_cur / cur_zoom + final_y;
            CxMin_cur = CxMin_cur / cur_zoom + final_x;
            CyMin_cur = CyMin_cur / cur_zoom + final_y;
        }

        /* */
        double PixelWidth = (CxMax_cur - CxMin_cur) / iXmax;
        double PixelHeight = (CyMax_cur - CyMin_cur) / iYmax;
        /* color component ( R or G or B) is coded from 0 to 255 */
        /* it is 24 bit color RGB file */
        const int MaxColorComponentValue = 255;
        unsigned int color[3];
        unsigned int ppmMatrix[iXmax][iYmax][3];
        for (iY = 0; iY < iYmax; iY++)
        {
            for (iX = 0; iX < iYmax; iX++)
            {
                for (int i = 0; i < 3; i++)
                {
                    ppmMatrix[iY][iX][i] = 0;
                }
            }
        }

        // Z=Zx+Zy*i  ;   Z0 = 0
        double Zx, Zy;
        double Zx2, Zy2; // Zx2=Zx*Zx;  Zy2=Zy*Zy
        //
        int Iteration;
        const int IterationMax = 2000;
        // bail-out value , radius of circle ;
        const double EscapeRadius = 300;
        double ER2 = EscapeRadius * EscapeRadius;

        //compute and write image data bytes to the file
        for (iY = frame_rank; iY < iYmax; iY += frame_size)
        {
            Cy = CyMin_cur + iY * PixelHeight;
            if (fabs(Cy) < PixelHeight / 2)
                Cy = 0.0; // Main antenna
            for (iX = 0; iX < iXmax; iX++)
            {
                Cx = CxMin_cur + iX * PixelWidth;
                // initial value of orbit = critical point Z= 0
                Zx = 0.0;
                Zy = 0.0;
                Zx2 = Zx * Zx;
                Zy2 = Zy * Zy;
                //
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

                for (int i = 0; i < 3; i++)
                {
                    ppmMatrix[iY][iX][i] = color[i];
                }
            }
        }

        //MPI_REDUCE
        unsigned int GLOBALppmMatrix[iXmax][iYmax][3];
        MPI_Reduce(ppmMatrix, GLOBALppmMatrix, iYmax * iXmax * 3, MPI_INT, MPI_SUM, frame_maser, frame_comm);

        //WRITE TO FILE
        if (my_rank == my_master)
        {
            FILE *fp;
            char filename[30];
            sprintf(filename, "zoom%d.ppm", actual_zoom);
            char *comment = "# "; //comment should start with #
            //create new file,give it a name and open it in binary mode
            fp = fopen(filename, "wb"); // b -  binary mode
            //write ASCII header to the file
            fprintf(fp, "P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);
            unsigned char CHARcolor[3];
            for (iY = 0; iY < iYmax; iY++)
            {
                for (iX = 0; iX < iYmax; iX++)
                {
                    //write color to the file
                    for (int i = 0; i < 3; i++)
                    {
                        CHARcolor[i] = GLOBALppmMatrix[iY][iX][i];
                    }
                    fwrite(CHARcolor, 1, 3, fp);
                }
            }

            fclose(fp);
        }

        // finalize the frame_comm comunicator
        MPI_Comm_free(&frame_comm);
    }

    MPI_Finalize();
    return 0;
}
