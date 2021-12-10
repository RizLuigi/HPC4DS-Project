#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
//#include <omp.h>
#include <math.h>
#include <argp.h>
#include <sys/stat.h>
#include <sys/types.h>

struct argp_option options[] =
    {
        {"threads", 't', "threadCount", 0, "Use threadCount threads per process"},
        {"zooms", 'z', "zoomCount", 0, "Create zoomCount images"},
        {"zoom-factor", 'f', "zoomFactor", 0, "Use zoomFactor as magnification factor"},
        {"final-x", 'x', "X", 0, "Use X as final x-coordinate"},
        {"final-y", 'y', "Y", 0, "Use Y as final y-coordinate"},
        {0}};

struct arguments
{
    int threads, zoom, zoom_factor;
    double final_x, final_y;
};

static int parse_opt(int key, char *arg, struct argp_state *state)
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

int mandelbrotIterations(double Cx, double Cy, int IterationMax, double ER2);
void colorFromIterations(int Iteration, int IterationMax, unsigned int *color);

int main(int argc, char *argv[])
{
    struct arguments arguments;

    // DEFAULTS
    arguments.threads = 10;
    arguments.zoom = 4;
    arguments.zoom_factor = 2;
    arguments.final_x = -1.1428;
    arguments.final_y = -0.3;

    argp_parse(&argp, argc, argv, 0, 0, &arguments);

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
    const int iXmax = 4000;
    const int iYmax = 4000;
    /* world ( double) coordinate = parameter plane*/
    double Cx, Cy;
    const double CxMin = -2.5;
    const double CxMax = 1.5;
    const double CyMin = -2.0;
    const double CyMax = 2.0;

    /* color component ( R or G or B) is coded from 0 to 255 */
    /* it is 24 bit color RGB file */
    const int MaxColorComponentValue = 255;

    //
    const int IterationMax = 2000;
    // bail-out value , radius of circle ;
    const double EscapeRadius = 300;
    double ER2 = EscapeRadius * EscapeRadius;

    if (comm_sz <= arguments.zoom)
    {
        unsigned int iterations[iXmax][iYmax];
        // In questo caso i processi sono indipendenti: ognuno lavora su uno o più zoom
        for (int i = 0; i < arguments.zoom / comm_sz; i++)
        {
            int actual_zoom = my_rank * arguments.zoom / comm_sz + i;
            printf("Process %d working on zoom %d\n", my_rank, actual_zoom);

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

            double PixelWidth = (CxMax_cur - CxMin_cur) / iXmax;
            double PixelHeight = (CyMax_cur - CyMin_cur) / iYmax;

            for (int iY = 0; iY < iYmax; iY++)
            {
                Cy = CyMin_cur + iY * PixelHeight;
                if (fabs(Cy) < PixelHeight / 2)
                    Cy = 0.0; // Main antenna
                for (int iX = 0; iX < iXmax; iX++)
                {
                    Cx = CxMin_cur + iX * PixelWidth;
                    iterations[iX][iY] = mandelbrotIterations(Cx, Cy, IterationMax, ER2);
                }
            }

            FILE *fp;
            mkdir("./frames", 0777);
            char filename[30];
            sprintf(filename, "frames/zoom%d.ppm", actual_zoom);
            char *comment = "# "; //comment should start with #
            //create new file,give it a name and open it in binary mode
            fp = fopen(filename, "wb"); // b -  binary mode
            //write ASCII header to the file
            fprintf(fp, "P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);
            unsigned int color[3];
            unsigned char CHARcolor[3];
            for (int iY = 0; iY < iYmax; iY++)
            {
                for (int iX = 0; iX < iXmax; iX++)
                {
                    colorFromIterations(iterations[iX][iY], IterationMax, color);
                    //write color to the file
                    for (int i = 0; i < 3; i++)
                    {
                        CHARcolor[i] = color[i];
                    }
                    fwrite(CHARcolor, 1, 3, fp);
                }
            }
        }
    }
    else
    {
        int proc_per_zoom = comm_sz / arguments.zoom;

        // creazione communicators gropus
        int working_frame = my_rank / proc_per_zoom;
        MPI_Comm frame_comm;
        MPI_Comm_split(MPI_COMM_WORLD, working_frame, comm_sz, &frame_comm);
        int frame_rank, frame_size;
        MPI_Comm_rank(frame_comm, &frame_rank);
        MPI_Comm_size(frame_comm, &frame_size);
        int frame_master = frame_rank - frame_rank % proc_per_zoom;

        //printf("%d/%d - frame_size: %d - frame_rank: %d frame_master: %d\n", my_rank, working_frame, frame_size, frame_rank, frame_master);

        // new coordinates for current frame
        double alpha = 1.0/pow(2.0, working_frame);
        double mean_x = alpha * ((CxMax + CxMin) / 2.0) + (1 - alpha) * final_x;
        double mean_y = alpha * ((CyMax + CyMin) / 2.0) + (1 - alpha) * final_y;
        double delta_x = ((CxMax - CxMin)/2.0) / pow(zoom_inc, working_frame);
        double delta_y = ((CyMax - CyMin)/2.0) / pow(zoom_inc, working_frame);

        double CxMax_cur = mean_x + delta_x;
        double CxMin_cur = mean_x - delta_x;
        double CyMax_cur = mean_y + delta_y;
        double CyMin_cur = mean_y - delta_y;

        printf("(%lf, %lf), (%lf, %lf)\n", mean_x, mean_y, delta_x, delta_y);

        /*if (working_frame != 0)
        {
            int cur_zoom = zoom_inc * working_frame;
            CxMax_cur = CxMax_cur / cur_zoom + final_x;
            CyMax_cur = CyMax_cur / cur_zoom + final_y;
            CxMin_cur = CxMin_cur / cur_zoom + final_x;
            CyMin_cur = CyMin_cur / cur_zoom + final_y;
        }*/

        /* */
        double PixelWidth = (CxMax_cur - CxMin_cur) / iXmax;
        double PixelHeight = (CyMax_cur - CyMin_cur) / iYmax;

        unsigned int rows_per_proc = iYmax / frame_size;
        unsigned int iterations[iXmax * rows_per_proc];

        for (int i = 0; i < iXmax * rows_per_proc; i++)
        {
            Cy = CyMin_cur + (i / iXmax + frame_rank * rows_per_proc) * PixelHeight;
            if (fabs(Cy) < PixelHeight / 2)
                Cy = 0.0; // Main antenna
            Cx = CxMin_cur + (i % iXmax) * PixelWidth;
            iterations[i] = mandelbrotIterations(Cx, Cy, IterationMax, ER2);
        }

        unsigned int *iterations_gathered;
        //MPI_GATHER
        if (frame_rank == frame_master)
        {
            iterations_gathered = (unsigned int *)malloc(iXmax * iYmax * sizeof(unsigned int));
            MPI_Gather(iterations, iXmax * rows_per_proc, MPI_INT, iterations_gathered, iXmax * rows_per_proc, MPI_INT, frame_master, frame_comm);
        }
        else
        {
            MPI_Gather(iterations, iXmax * rows_per_proc, MPI_INT, NULL, iXmax * rows_per_proc, MPI_INT, frame_master, frame_comm);
        }

        //WRITE TO FILE
        if (frame_rank == frame_master)
        {
            FILE *fp;
            mkdir("./frames", 0777);
            char filename[30];
            sprintf(filename, "frames/zoom%d.ppm", working_frame);
            char *comment = "# "; //comment should start with #
            //create new file,give it a name and open it in binary mode
            fp = fopen(filename, "wb"); // b -  binary mode
            //write ASCII header to the file
            fprintf(fp, "P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);
            unsigned int color[3];
            unsigned char CHARcolor[3];
            for (int i = 0; i < iXmax * iYmax; i++)
            {
                colorFromIterations(iterations_gathered[i], IterationMax, color);
                //write color to the file
                for (int i = 0; i < 3; i++)
                {
                    CHARcolor[i] = color[i];
                }
                fwrite(CHARcolor, 1, 3, fp);
            }

            fclose(fp);
        }

        // finalize the frame_comm comunicator
        MPI_Comm_free(&frame_comm);
    }

    MPI_Finalize();
    return 0;
}

int mandelbrotIterations(double Cx, double Cy, int IterationMax, double ER2)
{
    // initial value of orbit = critical point Z= 0
    double Zx = 0.0;
    double Zy = 0.0;
    double Zx2 = Zx * Zx;
    double Zy2 = Zy * Zy;
    //
    int Iteration;
    for (Iteration = 0; Iteration < IterationMax && ((Zx2 + Zy2) < ER2); Iteration++)
    {
        Zy = 2 * Zx * Zy + Cy;
        Zx = Zx2 - Zy2 + Cx;
        Zx2 = Zx * Zx;
        Zy2 = Zy * Zy;
    };
    return Iteration;
}

void colorFromIterations(int Iteration, int IterationMax, unsigned int *color)
{
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
    return;
}
