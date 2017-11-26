/* source: http://marathon.csee.usf.edu/edge/edge_detection.html */
/* URL: ftp://figment.csee.usf.edu/pub/Edge_Comparison/source_code/canny.src */

/* EECS 222 Assignment 5 solution */

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
typedef short int simg[SIZE];	/* define our communication data type */
typedef float float1[SIZE];	/* define our communication data type */
typedef float float2[WINSIZE];	/* define our communication data type */
typedef short int simg[SIZE];	/* define our communication data type */

DEFINE_I_TYPED_SENDER(img, img)		// creates interface i_img_sender
DEFINE_I_TYPED_RECEIVER(img, img)	// creates interface i_img_receiver
DEFINE_I_TYPED_TRANCEIVER(img, img)	// creates interface i_img_tranceiver

DEFINE_C_TYPED_QUEUE(img, img)		// creates channel c_img_queue


DEFINE_I_TYPED_SENDER(simg, simg)		// creates interface i_simg_sender
DEFINE_I_TYPED_RECEIVER(simg, simg)	// creates interface i_simg_receiver
DEFINE_I_TYPED_TRANCEIVER(simg, simg)	// creates interface i_simg_tranceiver

DEFINE_C_TYPED_QUEUE(simg, simg)		// creates channel c_simg_queue


DEFINE_I_TYPED_SENDER(float1, float1)		// creates interface i_float1_sender
DEFINE_I_TYPED_RECEIVER(float1, float1)	// creates interface i_float1_receiver
DEFINE_I_TYPED_TRANCEIVER(float1, float1)	// creates interface i_float1_tranceiver

DEFINE_C_TYPED_QUEUE(float1, float1)		// creates channel c_float1_queue

DEFINE_I_TYPED_SENDER(float2, float2)		// creates interface i_float2_sender
DEFINE_I_TYPED_RECEIVER(float2, float2)	// creates interface i_float2_receiver
DEFINE_I_TYPED_TRANCEIVER(float2, float2)	// creates interface i_float2_tranceiver

DEFINE_C_TYPED_QUEUE(float2, float2)		// creates channel c_float2_queue

behavior Stimulus(i_img_sender ImgOut)
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
	   }
	}
};

behavior Monitor(i_img_receiver ImgIn)
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

	   for(i=0; i<IMG_NUM; i++)
	   {
	      ImgIn.receive(&EdgeImage);

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

behavior receiveImage(i_img_receiver ImgIn, out unsigned char ImgOut[SIZE])
{
	unsigned char Image[SIZE];

	void main()
	{

	      ImgIn.receive(&Image);
	      ImgOut = Image;

	}
};

behavior make_gaussian_kernel(out float kernelout[WINSIZE]){
	
	int windowsize = WINSIZE;
	int center;
	float sigma = SIGMA;
	float kernel[WINSIZE];

	void make_gaussian_kernel1(float sigma, float *kernel, int *windowsize)
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

	void main(){

		if(VERBOSE) printf("   Computing the gaussian smoothing kernel.\n");
		make_gaussian_kernel1(sigma, kernel, &windowsize);

		//centerout = (*windowsize) / 2;
    kernelout = kernel;
	}

};

behavior blurX(in float kernelin[WINSIZE],in unsigned char ImgIn[SIZE],out float tempimout[SIZE])
{
	int r, c, rr, cc;
	int rows = ROWS,cols = COLS;
	float dot,            /* Dot product summing variable. */
	sum;
	float kernel[WINSIZE];
	unsigned char image[SIZE];
	int windowsize = WINSIZE;
	int center;
	float tempim[SIZE];

	void main()
	{

	kernel = kernelin;
  image = ImgIn;
	center = windowsize / 2;
		   /****************************************************************************
	   * Blur in the x - direction.
	   ****************************************************************************/
	   if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
	   for(r=0;r<rows;r++){
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

	tempimout = tempim;
	}// main ends
};

behavior blurY(in float kernelin[WINSIZE],in float tempimin[SIZE],out short int smoothedimout[SIZE])
{
	int r, c, rr, cc;
	int rows = ROWS,cols = COLS;
	float dot,            /* Dot product summing variable. */
	sum;
	float kernel[WINSIZE];
	short int smoothedim[SIZE];
	float tempim[SIZE];
	int windowsize = WINSIZE;
	int center;

	void main(){
	kernel = kernelin;
	tempim = tempimin;
	center = windowsize / 2;
		   /****************************************************************************
	   * Blur in the y - direction.
	   ****************************************************************************/
	   if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
	   for(c=0;c<cols;c++){
	      for(r=0;r<rows;r++){
	         sum = 0.0;
	         dot = 0.0;
	         for(rr=(-center);rr<=center;rr++){
	            if(((r+rr) >= 0) && ((r+rr) < rows)){
	               dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
	               sum += kernel[center+rr];
	            }
	         }
	         smoothedim[r*cols+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
	      }
	   }

	   smoothedimout = smoothedim;
		
	}//main ends
};



behavior gaussian_smooth(i_img_receiver ImgIn, out short sImgOut[SIZE])
{

  unsigned char Image[SIZE];
  float kernel[WINSIZE];
  float tempim[SIZE];

  receiveImage ri(ImgIn, Image);
	make_gaussian_kernel mgk(kernel);
 
	blurX bx (kernel,Image,tempim);
 
	blurY by (kernel,tempim,sImgOut);


	void main(){

  ri;
	mgk;
	bx;
	by;	

 }
};

behavior derivative_x_y(in short int sImgIn[SIZE], out short int xImgOut[SIZE],out short int yImgOut[SIZE]){
	short int smoothedim[SIZE];/* The input smooth image */
	short int delta_x[SIZE];    /* The output dx image */
	short int delta_y[SIZE];    /* The output dy image */
	int rows = ROWS;
	int cols = COLS;

	void derivative_x_y1(short int *smoothedim, int rows, int cols,
	        short int *delta_x, short int *delta_y)
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

	void main(){

		
		smoothedim = sImgIn;
		
		derivative_x_y1(smoothedim, rows, cols, delta_x, delta_y);

	   xImgOut = delta_x;
	   yImgOut = delta_y;
	   
	}
	};


behavior magnitude_x_y(in short int DxIn[SIZE],in short int DyIn[SIZE],out short int magOut[SIZE])
	{
		int r, c, pos, sq1, sq2;
		int rows = ROWS;
	    int cols = COLS;

		short int magnitude[SIZE];/* The output mag image */
		short int delta_x[SIZE];    /* The input dx image */
		short int delta_y[SIZE];    /* The input dy image */

		void main(){
		delta_x = DxIn;
		delta_y = DyIn;

		for(r=0,pos=0;r<rows;r++){
	      for(c=0;c<cols;c++,pos++){
	         sq1 = (int)delta_x[pos] * (int)delta_x[pos];
	         sq2 = (int)delta_y[pos] * (int)delta_y[pos];
	         magnitude[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
	      }
	   }

	   magOut = magnitude;
		}
	
	};


behavior non_max_supp(in short int magIn[SIZE],in short int DxIn[SIZE],in short int DyIn[SIZE],out unsigned char nmsOut[SIZE])
{
	int rowcount, colcount, count;
	short *magrowptr,*magptr;
	short *gxrowptr,*gxptr;
	short *gyrowptr,*gyptr,z1,z2;
	short m00,gx,gy;
	float mag1,mag2,xperp,yperp;
	unsigned char *resultrowptr, *resultptr;

	int rows = ROWS;
	int cols = COLS;
	short int magnitude[SIZE];/* The output mag image */
	short int delta_x[SIZE];    /* The input dx image */
	short int delta_y[SIZE];    /* The input dy image */
	unsigned char nms[SIZE];

	void non_max_supp1(short *mag, short *gradx, short *grady, int nrows, int ncols,
	    unsigned char *result)
	{
	    int rowcount, colcount,count;
	    short *magrowptr,*magptr;
	    short *gxrowptr,*gxptr;
	    short *gyrowptr,*gyptr,z1,z2;
	    short m00,gx,gy;
	    float mag1,mag2,xperp,yperp;
	    unsigned char *resultrowptr, *resultptr;

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

	void main(){

		delta_x = DxIn;
		delta_y = DyIn;
		magnitude = magIn;

		non_max_supp1(magnitude, delta_x, delta_y, rows, cols, nms);
		nmsOut = nms;
   
}  //Main ends
} ;  //behavior ends


behavior apply_hysteresis(in short int magIn[SIZE],in unsigned char nmsIn[SIZE],i_img_sender ImgOut)
{
	short int magnitude[SIZE];/* The output mag image */
	unsigned char nms[SIZE];
	unsigned char edge[SIZE];

	int rows = ROWS, cols = COLS;
	float tlow = TLOW, thigh = THIGH;

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
	void apply_hysteresis1(short int *mag, unsigned char *nms, int rows, int cols,
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

	void main(){
	
	
	magnitude = magIn;
	nms = nmsIn;

	apply_hysteresis1(magnitude, nms, rows, cols, tlow, thigh, edge);

	ImgOut.send(edge);


	}//main ends
};

behavior DUT(i_img_receiver ImgIn, i_img_sender ImgOut){

	unsigned char nmsq[SIZE];
	short int sml[SIZE],dx[SIZE],dy[SIZE],mag[SIZE];

	gaussian_smooth gs(ImgIn,sml);
	derivative_x_y dxy(sml,dx,dy);
	magnitude_x_y magxy(dx,dy,mag);
	non_max_supp nms(mag,dx,dy,nmsq);
	apply_hysteresis ahy(mag,nmsq,ImgOut);

	void main(){

	while(1){
	/****************************************************************************
	   * Perform gaussian smoothing on the image using the input standard
	   * deviation.
	   ****************************************************************************/
	   if(VERBOSE) printf("Smoothing the image using a gaussian kernel.\n");
	   gs;

	   /****************************************************************************
	   * Compute the first derivative in the x and y directions.
	   ****************************************************************************/
	   if(VERBOSE) printf("Computing the X and Y first derivatives.\n");
	   dxy;

	   /****************************************************************************
	   * Compute the magnitude of the gradient.
	   ****************************************************************************/
	   if(VERBOSE) printf("Computing the magnitude of the gradient.\n");
	   magxy;

	   /****************************************************************************
	   * Perform non-maximal suppression.
	   ****************************************************************************/
	   if(VERBOSE) printf("Doing the non-maximal suppression.\n");
	   nms;

	   /****************************************************************************
	   * Use hysteresis to mark the edge pixels.
	   ****************************************************************************/
	   if(VERBOSE) printf("Doing hysteresis thresholding.\n");
	   ahy;
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
	Stimulus stimulus(q1);
	Platform platform(q1, q2);
	Monitor monitor(q2);

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

