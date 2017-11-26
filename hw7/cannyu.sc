/* source: http://marathon.csee.usf.edu/edge/edge_detection.html */
/* URL: ftp://figment.csee.usf.edu/pub/Edge_Comparison/source_code/canny.src */

/* EECS 222 Assignment 6 solution */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <sim.sh>
#include <c_typed_queue.sh>	/* make the templates available */
#include <c_typed_double_handshake.sh>

#define VERBOSE 0

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0
#define BOOSTBLURFACTOR 90.0
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define SIGMA 0.6
#define TLOW  0.3
#define THIGH 0.8

#define COLS 2704
#define ROWS 1520
#define SIZE COLS*ROWS
#define VIDEONAME "EngPlaza"
#define IMG_IN    "video/" VIDEONAME "%03d.pgm"
#define IMG_OUT   VIDEONAME "%03d_edges.pgm"
#define IMG_NUM   20 /* number of images processed (1 or more) */
#define AVAIL_IMG 20 /* number of different image frames (1 or more) */


/* upper bound for the size of the gaussian kernel
 * SIGMA must be less than 4.0
 * check for 'windowsize' below
 */
#define WINSIZE 21

typedef unsigned char img[SIZE];	/* define our communication data type */
typedef short int     simg[SIZE];
typedef sim_time      time;
typedef sim_delta     del;
/*typedef sim_time_string   buffer;*/

DEFINE_I_TYPED_SENDER(img, img)		// creates interface i_img_sender
DEFINE_I_TYPED_RECEIVER(img, img)	// creates interface i_img_receiver
DEFINE_I_TYPED_TRANCEIVER(img, img)	// creates interface i_img_tranceiver
DEFINE_I_TYPED_SENDER(time, time)
DEFINE_I_TYPED_RECEIVER(time, time)
DEFINE_I_TYPED_TRANCEIVER(time, time)

DEFINE_C_TYPED_QUEUE(img, img)		// creates channel c_img_queue
DEFINE_C_TYPED_QUEUE(time, time)

behavior Stimulus(i_img_sender ImgOut, i_time_sender Timeout)
{
	unsigned char Image[SIZE];

	/******************************************************************************
	* Function: read_pgm_image
	* Purpose: This function reads in an image in PGM format. The image can be
	* read in from either a file or from standard input. The image is only read
	* from standard input when infilename = NULL. Because the PGM format includes
	* the number of columns and the number of rows in the image, these are read
	* from the file. Memory to store the image is allocated OUTSIDE this function.
	* The found image size is checked against the expected rows and cols.
	* All comments in the header are discarded in the process of reading the
	* image. Upon failure, this function returns 0, upon sucess it returns 1.
	******************************************************************************/
	int read_pgm_image(const char *infilename, unsigned char *image, int rows, int cols)
	{
	   FILE *fp;
	   char buf[71];
	   int r, c;

	   /***************************************************************************
	   * Open the input image file for reading if a filename was given. If no
	   * filename was provided, set fp to read from standard input.
	   ***************************************************************************/
	   if(infilename == NULL) fp = stdin;
	   else{
	      if((fp = fopen(infilename, "r")) == NULL){
	         fprintf(stderr, "Error reading the file %s in read_pgm_image().\n",
	            infilename);
	         return(0);
	      }
	   }

	   /***************************************************************************
	   * Verify that the image is in PGM format, read in the number of columns
	   * and rows in the image and scan past all of the header information.
	   ***************************************************************************/
	   fgets(buf, 70, fp);
	   if(strncmp(buf,"P5",2) != 0){
	      fprintf(stderr, "The file %s is not in PGM format in ", infilename);
	      fprintf(stderr, "read_pgm_image().\n");
	      if(fp != stdin) fclose(fp);
	      return(0);
	   }
	   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */
	   sscanf(buf, "%d %d", &c, &r);
	   if(c != cols || r != rows){
	      fprintf(stderr, "The file %s is not a %d by %d image in ", infilename, cols, rows);
	      fprintf(stderr, "read_pgm_image().\n");
	      if(fp != stdin) fclose(fp);
	      return(0);
	   }
	   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */

	   /***************************************************************************
	   * Read the image from the file.
	   ***************************************************************************/
	   if((unsigned)rows != fread(image, cols, rows, fp)){
	      fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
	      if(fp != stdin) fclose(fp);
	      return(0);
	   }

	   if(fp != stdin) fclose(fp);
	   return(1);
	}

	void main(void)
	{
	   int i=0, n=0;
	   char infilename[40];
	   sim_time t1;
	   sim_time_string  buffer1;
	   
	   for(i=0; i<IMG_NUM; i++)
	   {
	      n = i % AVAIL_IMG;
	      sprintf(infilename, IMG_IN, n+1);

	      /****************************************************************************
	      * Read in the image.
	      ****************************************************************************/
	      if(VERBOSE) printf("Reading the image %s.\n", infilename);
	      if(read_pgm_image(infilename, Image, ROWS, COLS) == 0){
	         fprintf(stderr, "Error reading the input image, %s.\n", infilename);
	         exit(1);
	      }
	      ImgOut.send(Image);
	      t1 = now();
	      printf("%s : Stimulus sent frame %d.\n",time2str(buffer1,t1/(1 MICRO_SEC)),i+1);
	      Timeout.send(t1);
	   }
	}
};

behavior Monitor(i_img_receiver ImgIn, i_time_receiver Timein)
{
	unsigned char EdgeImage[SIZE];

	/******************************************************************************
	* Function: write_pgm_image
	* Purpose: This function writes an image in PGM format. The file is either
	* written to the file specified by outfilename or to standard output if
	* outfilename = NULL. A comment can be written to the header if coment != NULL.
	******************************************************************************/
	int write_pgm_image(const char *outfilename, unsigned char *image, int rows,
	    int cols, const char *comment, int maxval)
	{
	   FILE *fp;

	   /***************************************************************************
	   * Open the output image file for writing if a filename was given. If no
	   * filename was provided, set fp to write to standard output.
	   ***************************************************************************/
	   if(outfilename == NULL) fp = stdout;
	   else{
	      if((fp = fopen(outfilename, "w")) == NULL){
	         fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
	            outfilename);
	         return(0);
	      }
	   }

	   /***************************************************************************
	   * Write the header information to the PGM file.
	   ***************************************************************************/
	   fprintf(fp, "P5\n%d %d\n", cols, rows);
	   if(comment != NULL)
	      if(strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
	   fprintf(fp, "%d\n", maxval);

	   /***************************************************************************
	   * Write the image data to the file.
	   ***************************************************************************/
	   if((unsigned)rows != fwrite(image, cols, rows, fp)){
	      fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
	      if(fp != stdout) fclose(fp);
	      return(0);
	   }

	   if(fp != stdout) fclose(fp);
	   return(1);
	}

	void main(void)
	{
	   char outfilename[128];    /* Name of the output "edge" image */
	   int i, n;
	   sim_time t2,t3,d;
	   sim_time_string buffer2;

	   for(i=0; i<IMG_NUM; i++)
	   {
	      ImgIn.receive(&EdgeImage);
	      Timein.receive(&t3);
	      t2 = now();
	      d = t2 - t3;
	      printf("%s : Monitor received frame",time2str(buffer2,t2/(1 MICRO_SEC)));
	      printf(" %d with %s us delay.\n",i+1,time2str(buffer2,d/(1 MICRO_SEC)));
	      
	      

	      /****************************************************************************
	      * Write out the edge image to a file.
	      ****************************************************************************/
	      n = i % AVAIL_IMG;
	      sprintf(outfilename, IMG_OUT, n+1);
	      if(VERBOSE) printf("Writing the edge image in the file %s.\n", outfilename);
	      if(write_pgm_image(outfilename, EdgeImage, ROWS, COLS,"", 255) == 0){
	         fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
	         exit(1);
	      }
	   }
	    printf("%s : Monitor exits simulation.\n",time2str(buffer2,t2));
	   if(VERBOSE) printf("Monitor exits simulation.\n");
	   sim_exit(0);	// done testing, quit the simulation
	}
};

behavior DataIn(i_img_receiver ImgIn, i_img_sender ImgOut)
{
	unsigned char Image[SIZE];

	void main()
	{
	   while(1)
	   {
	      ImgIn.receive(&Image);
	      ImgOut.send(Image);
	   }
	}
};

behavior DataOut(i_img_receiver ImgIn, i_img_sender ImgOut)
{
	unsigned char Image[SIZE];

	void main()
	{
	   while(1)
	   {
	      ImgIn.receive(&Image);
	      ImgOut.send(Image);
	   }
	}
};

behavior Receive_Image(i_img_receiver ImgIn, out img image)
{
	void main(void)
	{
	   img Image;
	   waitfor (0);
	   ImgIn.receive(&Image);
	   image = Image;
	}
};

behavior Gaussian_Kernel(out float gaussian_kernel[WINSIZE], out int kernel_center)
{
	/*******************************************************************************
	* PROCEDURE: make_gaussian_kernel
	* PURPOSE: Create a one dimensional gaussian kernel.
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/
	void make_gaussian_kernel(float sigma, float *kernel, int *windowsize)
	{
	   int i, center;
	   float x, fx, sum=0.0;

	   *windowsize = 1 + 2 * ceil(2.5 * sigma);
	   center = (*windowsize) / 2;

	   if(VERBOSE) printf("      The kernel has %d elements.\n", *windowsize);

	   for(i=0;i<(*windowsize);i++){
	      x = (float)(i - center);
	      fx = pow(2.71828, -0.5*x*x/(sigma*sigma)) / (sigma * sqrt(6.2831853));
	      kernel[i] = fx;
	      sum += fx;
	   }

	   for(i=0;i<(*windowsize);i++) kernel[i] /= sum;

	   if(VERBOSE){
	      printf("The filter coefficients are:\n");
	      for(i=0;i<(*windowsize);i++)
	         printf("kernel[%d] = %f\n", i, kernel[i]);
	   }
	}

	void main(void)
	{
	   int windowsize,       /* Dimension of the gaussian kernel. */
	       center;           /* Half of the windowsize. */
	   float kernel[WINSIZE] /* A one dimensional gaussian kernel. */
			= {0.0};
	   waitfor (0);

	   /****************************************************************************
	   * Create a 1-dimensional gaussian smoothing kernel.
	   ****************************************************************************/
	   if(VERBOSE) printf("   Computing the gaussian smoothing kernel.\n");
	   make_gaussian_kernel(SIGMA, kernel, &windowsize);
	   center = windowsize / 2;
	   gaussian_kernel = kernel;
	   kernel_center = center;
	}
};

behavior BlurX_Slice(in img image, in float kernel[WINSIZE], in int center, in int x, out float tempim[SIZE])
{
	void blur_x(int rows, int cols)
	{
	   int r, c, cc;         /* Counter variables. */
	   float dot,            /* Dot product summing variable. */
	         sum;            /* Sum of the kernel weights variable. */
	 /*  k1 = kernel;*/
	 /*  center1 = center;*/
	   /****************************************************************************
	   * Blur in the x - direction.
	   ****************************************************************************/
	   if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
	   for(r=(ROWS/8)*(rows-1);r<=(ROWS/8)*rows-1;r++){
	      for(c=0;c<cols;c++){
	         dot = 0.0;
	         sum = 0.0;
	         for(cc=(-center);cc<=center;cc++){
	            if(((c+cc) >= 0) && ((c+cc) < cols)){
	               dot += (float)image[r*cols+(c+cc)] * kernel[center+cc];
	               sum += kernel[center+cc];
	            }
	         }
	         tempim[r*cols+c] = dot/sum;
	      }
	   }
	}

	void main(void)
	{
	   waitfor (260.5 MILLI_SEC);
	   blur_x(x, COLS);
	}
};

behavior BlurX(in img image, in float kernel[WINSIZE], in int center,
                out float tempim[SIZE], out float k1[WINSIZE], out int center1)
{
	
	BlurX_Slice  sliceX1(image, kernel, center, 1, tempim);
	BlurX_Slice  sliceX2(image, kernel, center, 2, tempim);
	BlurX_Slice  sliceX3(image, kernel, center, 3, tempim);
	BlurX_Slice  sliceX4(image, kernel, center, 4, tempim);
	BlurX_Slice  sliceX5(image, kernel, center, 5, tempim);
	BlurX_Slice  sliceX6(image, kernel, center, 6, tempim);
	BlurX_Slice  sliceX7(image, kernel, center, 7, tempim);
	BlurX_Slice  sliceX8(image, kernel, center, 8, tempim);

	void main(void)
	{
		k1 = kernel;
		center1 = center;
		par{
		    sliceX1;
		    sliceX2;
		    sliceX3;
		    sliceX4;
		    sliceX5;
		    sliceX6;
		    sliceX7;
		    sliceX8;
	} 
	}			
};



behavior BlurY_Slice(in float tempim[SIZE], in float kernel[WINSIZE], in int center,in int y, out simg smoothedim)
{
	void blur_y(int rows, int cols)
	{
	   int r, c, rr;         /* Counter variables. */
	   float dot,            /* Dot product summing variable. */
	         sum;            /* Sum of the kernel weights variable. */

	   /****************************************************************************
	   * Blur in the y - direction.
	   ****************************************************************************/
	   if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
	   for(c=(COLS/8)*(cols-1);c<=(COLS/8)*cols-1;c++){
	      for(r=0;r<rows;r++){
	         sum = 0.0;
	         dot = 0.0;
	         for(rr=(-center);rr<=center;rr++){
	            if(((r+rr) >= 0) && ((r+rr) < rows)){
	               dot += tempim[(r+rr)*COLS+c] * kernel[center+rr];
	               sum += kernel[center+rr];
	            }
	         }
	         smoothedim[r*COLS+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
	      }
	   }
	}

	void main(void)
	{
	   waitfor (273.375 MILLI_SEC);
	   blur_y(ROWS, y);
	}
};


behavior BlurY(in float tempim[SIZE], in float k1[WINSIZE], in int center1, out simg smoothedim)
{
	BlurY_Slice  sliceY1(tempim, k1, center1, 1, smoothedim);
        BlurY_Slice  sliceY2(tempim, k1, center1, 2, smoothedim);
        BlurY_Slice  sliceY3(tempim, k1, center1, 3, smoothedim);
        BlurY_Slice  sliceY4(tempim, k1, center1, 4, smoothedim);
        BlurY_Slice  sliceY5(tempim, k1, center1, 5, smoothedim);
        BlurY_Slice  sliceY6(tempim, k1, center1, 6, smoothedim);
        BlurY_Slice  sliceY7(tempim, k1, center1, 7, smoothedim);
        BlurY_Slice  sliceY8(tempim, k1, center1, 8, smoothedim);

	void main(void)
        {       
                
                par{
                    sliceY1;
                    sliceY2;
                    sliceY3;
                    sliceY4;
                    sliceY5;
                    sliceY6;
                    sliceY7;
                    sliceY8;
	}
	}		
};

behavior Gaussian_Smooth(i_img_receiver ImgIn,out img image, out float kernel[WINSIZE], out int center)
{
	/*******************************************************************************
	* PROCEDURE: gaussian_smooth
	* PURPOSE: Blur an image with a gaussian filter.
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/

	Receive_Image   receive(ImgIn, image);
	Gaussian_Kernel gauss(kernel, center);
/*	BlurX           blurX(image, kernel, center, tempim);*/
/*	BlurY           blurY(tempim, kernel, center, smoothedim);*/

	void main(void)
	{
	   receive;
	   gauss;
	}
};

behavior Derivative_X_Y(in simg smoothedim, out simg delta_x, out simg delta_y)
{
	/*******************************************************************************
	* PROCEDURE: derivative_x_y
	* PURPOSE: Compute the first derivative of the image in both the x any y
	* directions. The differential filters that are used are:
	*
	*                                          -1
	*         dx =  -1 0 +1     and       dy =  0
	*                                          +1
	*
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/
	void derivative_x_y(int rows, int cols)
	{
	   int r, c, pos;

	   /****************************************************************************
	   * Compute the x-derivative. Adjust the derivative at the borders to avoid
	   * losing pixels.
	   ****************************************************************************/
	   if(VERBOSE) printf("   Computing the X-direction derivative.\n");
	   for(r=0;r<rows;r++){
	      pos = r * cols;
	      delta_x[pos] = smoothedim[pos+1] - smoothedim[pos];
	      pos++;
	      for(c=1;c<(cols-1);c++,pos++){
	         delta_x[pos] = smoothedim[pos+1] - smoothedim[pos-1];
	      }
	      delta_x[pos] = smoothedim[pos] - smoothedim[pos-1];
	   }

	   /****************************************************************************
	   * Compute the y-derivative. Adjust the derivative at the borders to avoid
	   * losing pixels.
	   ****************************************************************************/
	   if(VERBOSE) printf("   Computing the Y-direction derivative.\n");
	   for(c=0;c<cols;c++){
	      pos = c;
	      delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos];
	      pos += cols;
	      for(r=1;r<(rows-1);r++,pos+=cols){
	         delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos-cols];
	      }
	      delta_y[pos] = smoothedim[pos] - smoothedim[pos-cols];
	   }
	}

	void main(void)
	{
	   waitfor (267 MILLI_SEC);
	   derivative_x_y(ROWS, COLS);
	}
};

behavior Magnitude_X_Y(in simg delta_x, in simg delta_y, out simg magnitude, out simg delta_x1, out simg delta_y1)
{
	/*******************************************************************************
	* PROCEDURE: magnitude_x_y
	* PURPOSE: Compute the magnitude of the gradient. This is the square root of
	* the sum of the squared derivative values.
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/
	void magnitude_x_y(int rows, int cols)
	{
	   int r, c, pos, sq1, sq2;
	   delta_x1 = delta_x;
	   delta_y1 = delta_y;
	   for(r=0,pos=0;r<rows;r++){
	      for(c=0;c<cols;c++,pos++){
	         sq1 = (int)delta_x[pos] * (int)delta_x[pos];
	         sq2 = (int)delta_y[pos] * (int)delta_y[pos];
	         magnitude[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
	      }
	   }
	}

	void main(void)
	{
	   waitfor (267 MILLI_SEC);
	   magnitude_x_y(ROWS, COLS);
	}
};

behavior Non_Max_Supp(in simg gradx, in simg grady, in simg mag, out img nms, out simg magx)
{
	/*******************************************************************************
	* PROCEDURE: non_max_supp
	* PURPOSE: This routine applies non-maximal suppression to the magnitude of
	* the gradient image.
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/
	void non_max_supp(int nrows, int ncols, unsigned char *result)
	{
	    int rowcount, colcount,count;
	    short *magrowptr,*magptr;
	    short *gxrowptr,*gxptr;
	    short *gyrowptr,*gyptr,z1,z2;
	    short m00,gx,gy;
	    float mag1,mag2,xperp,yperp;
	    unsigned char *resultrowptr, *resultptr;
	    magx = mag;
	   /****************************************************************************
	   * Zero the edges of the result image.
	   ****************************************************************************/
	    for(count=0,resultrowptr=result,resultptr=result+ncols*(nrows-1);
	        count<ncols; resultptr++,resultrowptr++,count++){
	        *resultrowptr = *resultptr = (unsigned char) 0;
	    }

	    for(count=0,resultptr=result,resultrowptr=result+ncols-1;
	        count<nrows; count++,resultptr+=ncols,resultrowptr+=ncols){
	        *resultptr = *resultrowptr = (unsigned char) 0;
	    }

	   /****************************************************************************
	   * Suppress non-maximum points.
	   ****************************************************************************/
	   for(rowcount=1,magrowptr=mag+ncols+1,gxrowptr=gradx+ncols+1,
	      gyrowptr=grady+ncols+1,resultrowptr=result+ncols+1;
	      rowcount<=nrows-2;	// bug fix 3/29/17, RD
	      rowcount++,magrowptr+=ncols,gyrowptr+=ncols,gxrowptr+=ncols,
	      resultrowptr+=ncols){
	      for(colcount=1,magptr=magrowptr,gxptr=gxrowptr,gyptr=gyrowptr,
	         resultptr=resultrowptr;colcount<=ncols-2;	// bug fix 3/29/17, RD
	         colcount++,magptr++,gxptr++,gyptr++,resultptr++){
	         m00 = *magptr;
	         if(m00 == 0){
	            *resultptr = (unsigned char) NOEDGE;
	         }
	         else{
	            xperp = -(gx = *gxptr)/((float)m00);
	            yperp = (gy = *gyptr)/((float)m00);
	         }

	         if(gx >= 0){
	            if(gy >= 0){
	                    if (gx >= gy)
	                    {
	                        /* 111 */
	                        /* Left point */
	                        z1 = *(magptr - 1);
	                        z2 = *(magptr - ncols - 1);

	                        mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;

	                        /* Right point */
	                        z1 = *(magptr + 1);
	                        z2 = *(magptr + ncols + 1);

	                        mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
	                    }
	                    else
	                    {
	                        /* 110 */
	                        /* Left point */
	                        z1 = *(magptr - ncols);
	                        z2 = *(magptr - ncols - 1);

	                        mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

	                        /* Right point */
	                        z1 = *(magptr + ncols);
	                        z2 = *(magptr + ncols + 1);

	                        mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp;
	                    }
	                }
	                else
	                {
	                    if (gx >= -gy)
	                    {
	                        /* 101 */
	                        /* Left point */
	                        z1 = *(magptr - 1);
	                        z2 = *(magptr + ncols - 1);

	                        mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;

	                        /* Right point */
	                        z1 = *(magptr + 1);
	                        z2 = *(magptr - ncols + 1);

	                        mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
	                    }
	                    else
	                    {
	                        /* 100 */
	                        /* Left point */
	                        z1 = *(magptr + ncols);
	                        z2 = *(magptr + ncols - 1);

	                        mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

	                        /* Right point */
	                        z1 = *(magptr - ncols);
	                        z2 = *(magptr - ncols + 1);

	                        mag2 = (z1 - z2)*xperp  + (m00 - z1)*yperp;
	                    }
	                }
	            }
	            else
	            {
	                if ((gy = *gyptr) >= 0)
	                {
	                    if (-gx >= gy)
	                    {
	                        /* 011 */
	                        /* Left point */
	                        z1 = *(magptr + 1);
	                        z2 = *(magptr - ncols + 1);

	                        mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

	                        /* Right point */
	                        z1 = *(magptr - 1);
	                        z2 = *(magptr + ncols - 1);

	                        mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
	                    }
	                    else
	                    {
	                        /* 010 */
	                        /* Left point */
	                        z1 = *(magptr - ncols);
	                        z2 = *(magptr - ncols + 1);

	                        mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

	                        /* Right point */
	                        z1 = *(magptr + ncols);
	                        z2 = *(magptr + ncols - 1);

	                        mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
	                    }
	                }
	                else
	                {
	                    if (-gx > -gy)
	                    {
	                        /* 001 */
	                        /* Left point */
	                        z1 = *(magptr + 1);
	                        z2 = *(magptr + ncols + 1);

	                        mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

	                        /* Right point */
	                        z1 = *(magptr - 1);
	                        z2 = *(magptr - ncols - 1);

	                        mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
	                    }
	                    else
	                    {
	                        /* 000 */
	                        /* Left point */
	                        z1 = *(magptr + ncols);
	                        z2 = *(magptr + ncols + 1);

	                        mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

	                        /* Right point */
	                        z1 = *(magptr - ncols);
	                        z2 = *(magptr - ncols - 1);

	                        mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
	                    }
	                }
	            }

	            /* Now determine if the current point is a maximum point */

	            if ((mag1 > 0.0) || (mag2 > 0.0))
	            {
	                *resultptr = (unsigned char) NOEDGE;
	            }
	            else
	            {
	                if (mag2 == 0.0)
	                    *resultptr = (unsigned char) NOEDGE;
	                else
	                    *resultptr = (unsigned char) POSSIBLE_EDGE;
	            }
	        }
	    }
	}

	void main(void)
	{
	   
	   img result;
	   waitfor (1360 MILLI_SEC);
	   non_max_supp(ROWS, COLS, result);
	   nms = result;
	}
};

behavior Apply_Hysteresis(in simg mag, in img nms, i_img_sender ImgOut)
{
	unsigned char EdgeImage[SIZE];/* The output edge image */

	/*******************************************************************************
	* PROCEDURE: follow_edges
	* PURPOSE: This procedure edges is a recursive routine that traces edgs along
	* all paths whose magnitude values remain above some specifyable lower
	* threshhold.
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/
	void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval,
	   int cols)
	{
	   short *tempmagptr;
	   unsigned char *tempmapptr;
	   int i;
	   int x[8] = {1,1,0,-1,-1,-1,0,1},
	       y[8] = {0,1,1,1,0,-1,-1,-1};

	   for(i=0;i<8;i++){
	      tempmapptr = edgemapptr - y[i]*cols + x[i];
	      tempmagptr = edgemagptr - y[i]*cols + x[i];

	      if((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval)){
	         *tempmapptr = (unsigned char) EDGE;
	         follow_edges(tempmapptr,tempmagptr, lowval, cols);
	      }
	   }
	}

	/*******************************************************************************
	* PROCEDURE: apply_hysteresis
	* PURPOSE: This routine finds edges that are above some high threshhold or
	* are connected to a high pixel by a path of pixels greater than a low
	* threshold.
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/
	void apply_hysteresis(int rows, int cols,
		float tlow, float thigh, unsigned char *edge)
	{
	   int r, c, pos, numedges, highcount, lowthreshold, highthreshold, hist[32768];
	   short int maximum_mag;

	   /****************************************************************************
	   * Initialize the edge map to possible edges everywhere the non-maximal
	   * suppression suggested there could be an edge except for the border. At
	   * the border we say there can not be an edge because it makes the
	   * follow_edges algorithm more efficient to not worry about tracking an
	   * edge off the side of the image.
	   ****************************************************************************/
	   for(r=0,pos=0;r<rows;r++){
	      for(c=0;c<cols;c++,pos++){
		 if(nms[pos] == POSSIBLE_EDGE) edge[pos] = POSSIBLE_EDGE;
		 else edge[pos] = NOEDGE;
	      }
	   }

	   for(r=0,pos=0;r<rows;r++,pos+=cols){
	      edge[pos] = NOEDGE;
	      edge[pos+cols-1] = NOEDGE;
	   }
	   pos = (rows-1) * cols;
	   for(c=0;c<cols;c++,pos++){
	      edge[c] = NOEDGE;
	      edge[pos] = NOEDGE;
	   }

	   /****************************************************************************
	   * Compute the histogram of the magnitude image. Then use the histogram to
	   * compute hysteresis thresholds.
	   ****************************************************************************/
	   for(r=0;r<32768;r++) hist[r] = 0;
	   for(r=0,pos=0;r<rows;r++){
	      for(c=0;c<cols;c++,pos++){
		 if(edge[pos] == POSSIBLE_EDGE) hist[mag[pos]]++;
	      }
	   }

	   /****************************************************************************
	   * Compute the number of pixels that passed the nonmaximal suppression.
	   ****************************************************************************/
	   for(r=1,numedges=0;r<32768;r++){
	      if(hist[r] != 0) maximum_mag = r;
	      numedges += hist[r];
	   }

	   highcount = (int)(numedges * thigh + 0.5);

	   /****************************************************************************
	   * Compute the high threshold value as the (100 * thigh) percentage point
	   * in the magnitude of the gradient histogram of all the pixels that passes
	   * non-maximal suppression. Then calculate the low threshold as a fraction
	   * of the computed high threshold value. John Canny said in his paper
	   * "A Computational Approach to Edge Detection" that "The ratio of the
	   * high to low threshold in the implementation is in the range two or three
	   * to one." That means that in terms of this implementation, we should
	   * choose tlow ~= 0.5 or 0.33333.
	   ****************************************************************************/
	   r = 1;
	   numedges = hist[1];
	   while((r<(maximum_mag-1)) && (numedges < highcount)){
	      r++;
	      numedges += hist[r];
	   }
	   highthreshold = r;
	   lowthreshold = (int)(highthreshold * tlow + 0.5);

	   if(VERBOSE){
	      printf("The input low and high fractions of %f and %f computed to\n",
		 tlow, thigh);
	      printf("magnitude of the gradient threshold values of: %d %d\n",
		 lowthreshold, highthreshold);
	   }

	   /****************************************************************************
	   * This loop looks for pixels above the highthreshold to locate edges and
	   * then calls follow_edges to continue the edge.
	   ****************************************************************************/
	   for(r=0,pos=0;r<rows;r++){
	      for(c=0;c<cols;c++,pos++){
		 if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)){
	            edge[pos] = EDGE;
	            follow_edges((edge+pos), (mag+pos), lowthreshold, cols);
		 }
	      }
	   }

	   /****************************************************************************
	   * Set all the remaining possible edges to non-edges.
	   ****************************************************************************/
	   for(r=0,pos=0;r<rows;r++){
	      for(c=0;c<cols;c++,pos++) if(edge[pos] != EDGE) edge[pos] = NOEDGE;
	   }
	}

	void main(void)
	{
	   waitfor (382.5 MILLI_SEC);
	   apply_hysteresis(ROWS, COLS, TLOW, THIGH, EdgeImage);
	   ImgOut.send(EdgeImage);
	}
};

behavior DUT(i_img_receiver ImgIn, i_img_sender ImgOut)
{
	piped img image, nms;      /* Points that are local maximal magnitude. */
	piped simg smoothedim;   /* The image after gaussian smoothing.      */
	piped simg delta_x;      /* The first devivative image, x-direction. */
	piped simg delta_y;      /* The first derivative image, y-direction. */
	piped simg magnitude;    /* The magnitude of the gadient image.      */
	piped float kernel[WINSIZE], tempim[SIZE];
	piped piped float k1[WINSIZE];
	piped int   center;
 
 //**************************************
	piped piped int center1;
	piped piped simg delta_x1;
	piped piped simg delta_y1;
	piped piped simg magnitude1;
	
	
	Gaussian_Smooth  gaussian_smooth(ImgIn, image, kernel, center);
	BlurX            blurx(image,kernel, center, tempim, k1, center1);
	BlurY	         blury(tempim, k1, center1, smoothedim);  
	Derivative_X_Y   derivative_x_y(smoothedim, delta_x, delta_y);
	Magnitude_X_Y    magnitude_x_y(delta_x, delta_y, magnitude,delta_x1, delta_y1);
	Non_Max_Supp     non_max_supp(delta_x1, delta_y1, magnitude, nms, magnitude1);
	Apply_Hysteresis apply_hysteresis(magnitude1, nms, ImgOut);
		

	/*******************************************************************************
	* PROCEDURE: canny
	* PURPOSE: To perform canny edge detection.
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/
	void main(void)
	{
	   /****************************************************************************
	   * Perform the edge detection. All of the work takes place here.
	   ****************************************************************************/
	  
           int i;	
	   pipe(i=0;i<20;i++)
	   {
	      
	      /****************************************************************************
	      * Perform gaussian smoothing on the image using the input standard
	      * deviation.
	      ****************************************************************************/
	      gaussian_smooth;

	      /****************************************************************************
	      * Compute the first derivative in the x and y directions.
	      ****************************************************************************/
	

	      blurx;

	      blury;		

	      derivative_x_y;

	      /****************************************************************************
	      * Compute the magnitude of the gradient.
	      ****************************************************************************/
	      magnitude_x_y;

	      /****************************************************************************
	      * Perform non-maximal suppression.
	      ****************************************************************************/
	      non_max_supp;

	      /****************************************************************************
	      * Use hysteresis to mark the edge pixels.
	      ****************************************************************************/
	      apply_hysteresis;
		
	   }
	
	}
};

behavior Platform(i_img_receiver ImgIn, i_img_sender ImgOut)
{
	c_img_queue q1(2ul), q2(2ul);
	DataIn din(ImgIn, q1);
	DUT canny(q1, q2);
	DataOut dout(q2, ImgOut);

	void main()
	{
	   par{
	      din;
	      canny;
	      dout;
	   }
	}
};


behavior Main(void)
{
	c_img_queue q1(2ul), q2(2ul);
        c_time_queue q3(5ul);
	Stimulus stimulus(q1,q3);
	Platform platform(q1, q2);
	Monitor monitor(q2,q3);

	int main(void)
	{
	   par{
	      stimulus;
	      platform;
	      monitor;
	   }
	   return 0; // never reached
	}
};

