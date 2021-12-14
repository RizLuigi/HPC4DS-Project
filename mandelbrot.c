#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <argp.h>
#include <sys/stat.h>
#include <sys/types.h>
 #include <unistd.h>

struct argp_option options[] =
    {
        {"threads", 't', "threadCount", 0, "Use threadCount threads per process"},
        {"zooms", 'z', "zoomCount", 0, "Create zoomCount images"},
        {"zoom-factor", 'f', "zoomFactor", 0, "Use zoomFactor as magnification factor"},
        {"final-x", 'x', "X", 0, "Use X as final x-coordinate"},
        {"final-y", 'y', "Y", 0, "Use Y as final y-coordinate"},
        {"resolution", 'r', "R", 0, "Create a R x R image"},
        {"iterations", 'i', "ITER", 0, "Use ITER iterations when evaluating the membership of a point to the Mandelbrot set"},
        {0}};

struct arguments
{
    int threads, zoom, zoom_factor, resolution, iterations;
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
    case 'r':
        arguments->resolution = atoi(arg);
        break;
    case 'i':
        arguments->iterations = atoi(arg);
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
    arguments.final_x = -1.768778833;
    arguments.final_y = -0.001738996;
    arguments.resolution = 400;
    arguments.resolution = 2000;

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
    const int iXmax = arguments.resolution;
    const int iYmax = arguments.resolution;
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
    const int IterationMax = arguments.iterations;
    // bail-out value , radius of circle ;
    const double EscapeRadius = 2.0;
    double ER2 = EscapeRadius * EscapeRadius;

    double start, finish, elapsed;
    start = omp_get_wtime();

    char foldername[50];
    if (my_rank == 0) {
        //srand(time(NULL));
        int r = getpid();
        sprintf(foldername, "./frames%d", r);
        mkdir(foldername, 0777);
    }

    //BROADCAST foldername
    MPI_Bcast( foldername , 50 , MPI_CHAR , 0 , MPI_COMM_WORLD);
    

    if (comm_sz <= arguments.zoom)
    {   
        //Same as: unsigned int iterations[iXmax][iYmax];
        unsigned int (*iterations)[iYmax];
        iterations = malloc(sizeof(int[iXmax][iYmax]));

        // In questo caso i processi sono indipendenti: ognuno lavora su uno o più zoom
        for (int i = 0; i < arguments.zoom / comm_sz; i++)
        {
            int actual_zoom = my_rank * arguments.zoom / comm_sz + i;
            //printf("Process %d working on zoom %d\n", my_rank, actual_zoom);

            // new coordinates for current frame
            double alpha = 1.0 / pow(2.0, actual_zoom);
            double mean_x = alpha * ((CxMax + CxMin) / 2.0) + (1 - alpha) * final_x;
            double mean_y = alpha * ((CyMax + CyMin) / 2.0) + (1 - alpha) * final_y;
            double delta_x = ((CxMax - CxMin) / 2.0) / pow(zoom_inc, actual_zoom);
            double delta_y = ((CyMax - CyMin) / 2.0) / pow(zoom_inc, actual_zoom);

            double CxMax_cur = mean_x + delta_x;
            double CxMin_cur = mean_x - delta_x;
            double CyMax_cur = mean_y + delta_y;
            double CyMin_cur = mean_y - delta_y;

            //printf("%d (%lf, %lf), (%lf, %lf)\n", actual_zoom, mean_x, mean_y, delta_x, delta_y);

            double PixelWidth = (CxMax_cur - CxMin_cur) / iXmax;
            double PixelHeight = (CyMax_cur - CyMin_cur) / iYmax;
            #ifdef _OPENMP
            #pragma omp parallel for num_threads(arguments.threads) private(Cy, Cx)
            #endif
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
            char filename[100];
            sprintf(filename, "%s/zoom%d.ppm", foldername, actual_zoom);
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

        // new coordinates for current frame
        double alpha = 1.0 / pow(2.0, working_frame);
        double mean_x = alpha * ((CxMax + CxMin) / 2.0) + (1 - alpha) * final_x;
        double mean_y = alpha * ((CyMax + CyMin) / 2.0) + (1 - alpha) * final_y;
        double delta_x = ((CxMax - CxMin) / 2.0) / pow(zoom_inc, working_frame);
        double delta_y = ((CyMax - CyMin) / 2.0) / pow(zoom_inc, working_frame);

        double CxMax_cur = mean_x + delta_x;
        double CxMin_cur = mean_x - delta_x;
        double CyMax_cur = mean_y + delta_y;
        double CyMin_cur = mean_y - delta_y;

        //printf("(%lf, %lf), (%lf, %lf)\n", mean_x, mean_y, delta_x, delta_y);

        /* */
        double PixelWidth = (CxMax_cur - CxMin_cur) / iXmax;
        double PixelHeight = (CyMax_cur - CyMin_cur) / iYmax;

        unsigned int rows_per_proc = iYmax / frame_size;

        unsigned int (*iterations);
        iterations = malloc(sizeof(int[iXmax * rows_per_proc]));
        //unsigned int iterations[iXmax * rows_per_proc];

        #ifdef _OPENMP
        #pragma omp parallel for num_threads(arguments.threads) private(Cy, Cx)
        #endif
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
            //riempire le righe mancanti
            unsigned int rem_rows = iYmax % frame_size;
            //printf("Righe in eccesso: %d\n", rem_rows);
            for (int i = 0; i < iXmax * rem_rows; i++)
            {
                Cy = CyMin_cur + ((iYmax - rem_rows) + i / iXmax) * PixelHeight;
                if (fabs(Cy) < PixelHeight / 2)
                    Cy = 0.0; // Main antenna
                Cx = CxMin_cur + (i % iXmax) * PixelWidth;
                iterations_gathered[(iYmax - rem_rows) * iXmax + i] = mandelbrotIterations(Cx, Cy, IterationMax, ER2);
            }
        }
        else
        {
            MPI_Gather(iterations, iXmax * rows_per_proc, MPI_INT, NULL, iXmax * rows_per_proc, MPI_INT, frame_master, frame_comm);
        }

        free(iterations);

        //WRITE TO FILE
        if (frame_rank == frame_master)
        {
            FILE *fp;
            char filename[30];
            sprintf(filename, "%s/zoom%d.ppm", foldername, working_frame);
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

    MPI_Barrier(MPI_COMM_WORLD);

    if (my_rank==0) {
        finish = omp_get_wtime();
        elapsed = finish - start;
        printf("Elapsed time = %e seconds\n", elapsed);
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
