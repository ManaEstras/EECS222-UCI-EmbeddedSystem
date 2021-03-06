//////////////////////////////////////////////////////////////////////
// C++ source file generated by SpecC V2.2.1
// Design: Canny5
// File:   Canny5.cc
// Time:   Sat May 20 22:04:49 2017
//////////////////////////////////////////////////////////////////////

// Note: User-defined include files are inlined in this file.

// Note: System-defined include files are inlined in this file.

#include "Canny5.h"


unsigned int _IDcnt = 0;
// channel class definitions /////////////////////////////////////////

c_img_queue::c_img_queue(const unsigned long int (&size))
    : _specc::channel(), size(size),
    buffer(0),
    n(0ul),
    p(0ul),
    wr(0ul),
    ws(0ul)
{   
}

c_img_queue::~c_img_queue(void)
{   
}

#line 49 "Canny5.sc"
void c_img_queue::cleanup(void) { if ( !n) { free(buffer); buffer = 0;
    }
}

#line 49 "Canny5.sc"
void c_img_queue::receive(unsigned char (*d)[4110080]) { while( !n) { wr++ ; _specc::wait(event(&r), ((void*)0)); wr-- ;
    }

#line 49 "Canny5.sc"
    if (n <= p) { { unsigned int _scc_index_0; for(_scc_index_0=0;_scc_index_0<4110080;_scc_index_0++) ( *d)[_scc_index_0] = (buffer[p - n])[_scc_index_0]; }
    }
    else 

#line 49 "Canny5.sc"
    {    { unsigned int _scc_index_0; for(_scc_index_0=0;_scc_index_0<4110080;_scc_index_0++) ( *d)[_scc_index_0] = (buffer[p + size - n])[_scc_index_0]; }
    }

#line 49 "Canny5.sc"
    n-- ; if (ws) { _specc::notify(event(&s), ((void*)0));
    }

#line 49 "Canny5.sc"
    cleanup();
}

#line 49 "Canny5.sc"
void c_img_queue::send(unsigned char d[4110080]) { while(n >= size) { ws++ ; _specc::wait(event(&s), ((void*)0)); ws-- ;
    }

#line 49 "Canny5.sc"
    setup(); { unsigned int _scc_index_0; for(_scc_index_0=0;_scc_index_0<4110080;_scc_index_0++) (buffer[p])[_scc_index_0] = (d)[_scc_index_0]; } p++ ; if (p >= size) { p = 0;
    }

#line 49 "Canny5.sc"
    n++ ; if (wr) { _specc::notify(event(&r), ((void*)0));
    }
}

#line 49 "Canny5.sc"
void c_img_queue::setup(void) { if ( !buffer) { unsigned char dummy[4110080]; unsigned long int i; if ( !(buffer = (unsigned char (*)[4110080])malloc(sizeof(unsigned char [4110080]) * size))) { perror("c_typed_queue"); abort();
	}

#line 49 "Canny5.sc"
	for(i = 0; i < size; i++ ) { memcpy( &buffer[i],  &dummy, sizeof(unsigned char [4110080]));
	}
    }
}

// behavior class definitions ////////////////////////////////////////

#line 84 "Canny5.cc"
Stimulus::Stimulus(unsigned int _idcnt, i_img_sender (&ImgOut))
    : _specc::behavior(_idcnt), ImgOut(ImgOut)
{   
}

Stimulus::~Stimulus(void)
{   
}

#line 119 "Canny5.sc"
void Stimulus::main(void)
{   
    int i = 0; int n = 0;
    char infilename[40];

    for(i = 0; i < 20; i++ )
    {   
	n = i % 20;
	sprintf(infilename, "video/EngPlaza%03d.pgm", n + 1);




	if (0) printf("Reading the image %s.\n", infilename);
	if (read_pgm_image(infilename, Image, 1520, 2704) == 0) {
	    fprintf(stderr, "Error reading the input image, %s.\n", infilename);
	    exit(1);
	}
	ImgOut.send(Image);
    }
}

#line 66 "Canny5.sc"
int Stimulus::read_pgm_image(const char *infilename, unsigned char *image, int rows, int cols)
{   
    struct _IO_FILE *fp;
    char buf[71];
    int c; int r;

#line 76 "Canny5.sc"
    if (infilename == ((void *)0)) fp = stdin;
    else  {
	if ((fp = fopen(infilename, "r")) == ((void *)0)) {
	    fprintf(stderr, "Error reading the file %s in read_pgm_image().\n", 
		infilename);
	    return (0);
	}
    }

#line 89 "Canny5.sc"
    fgets(buf, 70, fp);
    if (strncmp(buf, "P5", 2) != 0) {
	fprintf(stderr, "The file %s is not in PGM format in ", infilename);
	fprintf(stderr, "read_pgm_image().\n");
	if (fp != stdin) fclose(fp);
	return (0);
    }
    do  { fgets(buf, 70, fp);
    }
    while(

#line 96 "Canny5.sc"
	buf[0] == '#');
    __isoc99_sscanf(buf, "%d %d",  &c,  &r);
    if (c != cols || r != rows) {
	fprintf(stderr, "The file %s is not a %d by %d image in ", infilename, cols, rows);
	fprintf(stderr, "read_pgm_image().\n");
	if (fp != stdin) fclose(fp);
	return (0);
    }
    do  { fgets(buf, 70, fp);
    }
    while(

#line 104 "Canny5.sc"
	buf[0] == '#');




    if ((unsigned int)rows != fread(image, cols, rows, fp)) {
	fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
	if (fp != stdin) fclose(fp);
	return (0);
    }

    if (fp != stdin) fclose(fp);
    return (1);
}

#line 175 "Canny5.cc"
Monitor::Monitor(unsigned int _idcnt, i_img_receiver (&ImgIn))
    : _specc::behavior(_idcnt), ImgIn(ImgIn)
{   
}

Monitor::~Monitor(void)
{   
}

#line 191 "Canny5.sc"
void Monitor::main(void)
{   
    char outfilename[128];
    int i; int n;

    for(i = 0; i < 20; i++ )
    {   
	ImgIn.receive( &EdgeImage);




	n = i % 20;
	sprintf(outfilename, "EngPlaza%03d_edges.pgm", n + 1);
	if (0) printf("Writing the edge image in the file %s.\n", outfilename);
	if (write_pgm_image(outfilename, EdgeImage, 1520, 2704, "", 255) == 0) {
	    fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
	    exit(1);
	}
    }
    if (0) printf("Monitor exits simulation.\n");
    sim_exit(0);
}

#line 152 "Canny5.sc"
int Monitor::write_pgm_image(const char *outfilename, unsigned char *image, int rows, 
    int cols, const char *comment, int maxval)
{   
    struct _IO_FILE *fp;

#line 161 "Canny5.sc"
    if (outfilename == ((void *)0)) fp = stdout;
    else  {
	if ((fp = fopen(outfilename, "w")) == ((void *)0)) {
	    fprintf(stderr, "Error writing the file %s in write_pgm_image().\n", 
		outfilename);
	    return (0);
	}
    }




    fprintf(fp, "P5\n%d %d\n", cols, rows);
    if (comment != ((void *)0))
	if (strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
    fprintf(fp, "%d\n", maxval);




    if ((unsigned int)rows != fwrite(image, cols, rows, fp)) {
	fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
	if (fp != stdout) fclose(fp);
	return (0);
    }

    if (fp != stdout) fclose(fp);
    return (1);
}

#line 247 "Canny5.cc"
DataIn::DataIn(unsigned int _idcnt, i_img_receiver (&ImgIn), i_img_sender (&ImgOut))
    : _specc::behavior(_idcnt), ImgIn(ImgIn), ImgOut(ImgOut)
{   
}

DataIn::~DataIn(void)
{   
}

#line 220 "Canny5.sc"
void DataIn::main()
{   
    while(1)
    {   
	ImgIn.receive( &Image);
	ImgOut.send(Image);
    }
}

#line 267 "Canny5.cc"
DataOut::DataOut(unsigned int _idcnt, i_img_receiver (&ImgIn), i_img_sender (&ImgOut))
    : _specc::behavior(_idcnt), ImgIn(ImgIn), ImgOut(ImgOut)
{   
}

DataOut::~DataOut(void)
{   
}

#line 234 "Canny5.sc"
void DataOut::main()
{   
    while(1)
    {   
	ImgIn.receive( &Image);
	ImgOut.send(Image);
    }
}

#line 287 "Canny5.cc"
DUT::DUT(unsigned int _idcnt, i_img_receiver (&ImgIn), i_img_sender (&ImgOut))
    : _specc::behavior(_idcnt), ImgIn(ImgIn), ImgOut(ImgOut)
{   
}

DUT::~DUT(void)
{   
}

#line 655 "Canny5.sc"
void DUT::apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols, 
    float tlow, float thigh, unsigned char *edge)
{   
    int c; int highcount; int highthreshold; int hist[32768]; int lowthreshold; int numedges; int pos; int r;
    short int maximum_mag;

#line 668 "Canny5.sc"
    for(r = 0 , pos = 0; r < rows; r++ ) {
	for(c = 0; c < cols; c++  , pos++ ) {
	    if (nms[pos] == 128) edge[pos] = 128;
	    else  edge[pos] = 255;
	}
    }

    for(r = 0 , pos = 0; r < rows; r++  , pos += cols) {
	edge[pos] = 255;
	edge[pos + cols - 1] = 255;
    }
    pos = (rows - 1) * cols;
    for(c = 0; c < cols; c++  , pos++ ) {
	edge[c] = 255;
	edge[pos] = 255;
    }

#line 689 "Canny5.sc"
    for(r = 0; r < 32768; r++ ) hist[r] = 0;
    for(r = 0 , pos = 0; r < rows; r++ ) {
	for(c = 0; c < cols; c++  , pos++ ) {
	    if (edge[pos] == 128) hist[mag[pos]]++ ;
	}
    }




    for(r = 1 , numedges = 0; r < 32768; r++ ) {
	if (hist[r] != 0) maximum_mag = r;
	numedges += hist[r];
    }

    highcount = (int)(numedges * thigh + 5.000000000000000e-01);

#line 716 "Canny5.sc"
    r = 1;
    numedges = hist[1];
    while((r < (maximum_mag - 1)) && (numedges < highcount)) {
	r++ ;
	numedges += hist[r];
    }
    highthreshold = r;
    lowthreshold = (int)(highthreshold * tlow + 5.000000000000000e-01);

    if (0) {
	printf("The input low and high fractions of %f and %f computed to\n", 
	    tlow, thigh);
	printf("magnitude of the gradient threshold values of: %d %d\n", 
	    lowthreshold, highthreshold);
    }

#line 736 "Canny5.sc"
    for(r = 0 , pos = 0; r < rows; r++ ) {
	for(c = 0; c < cols; c++  , pos++ ) {
	    if ((edge[pos] == 128) && (mag[pos] >= highthreshold)) {
		edge[pos] = 0;
		follow_edges((edge + pos), (mag + pos), lowthreshold, cols);
	    }
	}
    }




    for(r = 0 , pos = 0; r < rows; r++ ) {
	for(c = 0; c < cols; c++  , pos++ ) if (edge[pos] != 0) edge[pos] = 255;
    }
}

#line 759 "Canny5.sc"
void DUT::canny(unsigned char *image, int rows, int cols, float sigma, 
    float tlow, float thigh, unsigned char *edge)
{   
    unsigned char nms[4110080] = { 
      ((unsigned char)'\000') };
    short int smoothedim[4110080] = { 
      ((short int)0) };
    short int delta_x[4110080] = { 
      ((short int)0) };
    short int delta_y[4110080] = { 
      ((short int)0) };
    short int magnitude[4110080] = { 
      ((short int)0) };

#line 777 "Canny5.sc"
    if (0) printf("Smoothing the image using a gaussian kernel.\n");
    gaussian_smooth(image, rows, cols, sigma, smoothedim);




    if (0) printf("Computing the X and Y first derivatives.\n");
    derivative_x_y(smoothedim, rows, cols, delta_x, delta_y);




    if (0) printf("Computing the magnitude of the gradient.\n");
    magnitude_x_y(delta_x, delta_y, rows, cols, magnitude);




    if (0) printf("Doing the non-maximal suppression.\n");
    non_max_supp(magnitude, delta_x, delta_y, rows, cols, nms);




    if (0) printf("Doing hysteresis thresholding.\n");
    apply_hysteresis(magnitude, nms, rows, cols, tlow, thigh, edge);
}

#line 356 "Canny5.sc"
void DUT::derivative_x_y(short int *smoothedim, int rows, int cols, 
    short int *delta_x, short int *delta_y)
{   
    int c; int pos; int r;

#line 365 "Canny5.sc"
    if (0) printf("   Computing the X-direction derivative.\n");
    for(r = 0; r < rows; r++ ) {
	pos = r * cols;
	delta_x[pos] = smoothedim[pos + 1] - smoothedim[pos];
	pos++ ;
	for(c = 1; c < (cols - 1); c++  , pos++ ) {
	    delta_x[pos] = smoothedim[pos + 1] - smoothedim[pos - 1];
	}
	delta_x[pos] = smoothedim[pos] - smoothedim[pos - 1];
    }

#line 380 "Canny5.sc"
    if (0) printf("   Computing the Y-direction derivative.\n");
    for(c = 0; c < cols; c++ ) {
	pos = c;
	delta_y[pos] = smoothedim[pos + cols] - smoothedim[pos];
	pos += cols;
	for(r = 1; r < (rows - 1); r++  , pos += cols) {
	    delta_y[pos] = smoothedim[pos + cols] - smoothedim[pos - cols];
	}
	delta_y[pos] = smoothedim[pos] - smoothedim[pos - cols];
    }
}

#line 627 "Canny5.sc"
void DUT::follow_edges(unsigned char *edgemapptr, short int *edgemagptr, short int lowval, 
    int cols)
{   
    short int *tempmagptr;
    unsigned char *tempmapptr;
    int i;
    int x[8] = { 1,1,0,-1,-1,-1,0,1 };
    int y[8] = { 0,1,1,1,0,-1,-1,-1 };

    for(i = 0; i < 8; i++ ) {
	tempmapptr = edgemapptr - y[i] * cols + x[i];
	tempmagptr = edgemagptr - y[i] * cols + x[i];

	if (( *tempmapptr == 128) && ( *tempmagptr > lowval)) {
	     *tempmapptr = (unsigned char)0;
	    follow_edges(tempmapptr, tempmagptr, lowval, cols);
	}
    }
}

#line 287 "Canny5.sc"
void DUT::gaussian_smooth(unsigned char *image, int rows, int cols, float sigma, 
    short int *smoothedim)
{   
    int c; int cc; int r; int rr;
    int windowsize;
    int center;
    float tempim[4110080] = { 
      0.000000e+00f };
    float kernel[21] = { 
      0.000000e+00f };
    float dot;
    float sum;




    if (0) printf("   Computing the gaussian smoothing kernel.\n");
    make_gaussian_kernel(sigma, kernel,  &windowsize);
    center = windowsize / 2;




    if (0) printf("   Bluring the image in the X-direction.\n");
    for(r = 0; r < rows; r++ ) {
	for(c = 0; c < cols; c++ ) {
	    dot = 0.000000000000000e+00;
	    sum = 0.000000000000000e+00;
	    for(cc = ( -center); cc <= center; cc++ ) {
		if (((c + cc) >= 0) && ((c + cc) < cols)) {
		    dot += (float)image[r * cols + (c + cc)] * kernel[center + cc];
		    sum += kernel[center + cc];
		}
	    }
	    tempim[r * cols + c] = dot / sum;
	}
    }




    if (0) printf("   Bluring the image in the Y-direction.\n");
    for(c = 0; c < cols; c++ ) {
	for(r = 0; r < rows; r++ ) {
	    sum = 0.000000000000000e+00;
	    dot = 0.000000000000000e+00;
	    for(rr = ( -center); rr <= center; rr++ ) {
		if (((r + rr) >= 0) && ((r + rr) < rows)) {
		    dot += tempim[(r + rr) * cols + c] * kernel[center + rr];
		    sum += kernel[center + rr];
		}
	    }
	    smoothedim[r * cols + c] = (short int)(dot * 9.000000000000000e+01 / sum + 5.000000000000000e-01);
	}
    }
}

#line 399 "Canny5.sc"
void DUT::magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols, 
    short int *magnitude)
{   
    int c; int pos; int r; int sq1; int sq2;

    for(r = 0 , pos = 0; r < rows; r++ ) {
	for(c = 0; c < cols; c++  , pos++ ) {
	    sq1 = (int)delta_x[pos] * (int)delta_x[pos];
	    sq2 = (int)delta_y[pos] * (int)delta_y[pos];
	    magnitude[pos] = (short int)(5.000000000000000e-01 + sqrt((float)sq1 + (float)sq2));
	}
    }
}

#line 805 "Canny5.sc"
void DUT::main(void)
{   
    while(1)
    {   
	ImgIn.receive( &Image);




	if (0) printf("Starting Canny edge detection.\n");
	canny(Image, 1520, 2704, 6.000000000000000e-01, 3.000000000000000e-01, 8.000000000000000e-01, EdgeImage);

	ImgOut.send(EdgeImage);
    }
}

#line 255 "Canny5.sc"
void DUT::make_gaussian_kernel(float sigma, float *kernel, int *windowsize)
{   
    int center; int i;
    float fx; float sum = 0.000000e+00f; float x;

     *windowsize = 1 + 2 * ceil(2.500000000000000e+00 * sigma);
    center = ( *windowsize) / 2;

    if (0) printf("      The kernel has %d elements.\n",  *windowsize);

    for(i = 0; i < ( *windowsize); i++ ) {
	x = (float)(i - center);
	fx = pow(2.718280000000000e+00,  -5.000000000000000e-01 * x * x / (sigma * sigma)) / (sigma * sqrt(6.283185300000000e+00));
	kernel[i] = fx;
	sum += fx;
    }

    for(i = 0; i < ( *windowsize); i++ ) kernel[i] /= sum;

    if (0) {
	printf("The filter coefficients are:\n");
	for(i = 0; i < ( *windowsize); i++ )
	    printf("kernel[%d] = %f\n", i, kernel[i]);
    }
}

#line 421 "Canny5.sc"
void DUT::non_max_supp(short int *mag, short int *gradx, short int *grady, int nrows, int ncols, 
    unsigned char *result)
{   
    int colcount; int count; int rowcount;
    short int *magptr; short int *magrowptr;
    short int *gxptr; short int *gxrowptr;
    short int *gyptr; short int *gyrowptr; short int z1; short int z2;
    short int gx; short int gy; short int m00;
    float mag1; float mag2; float xperp; float yperp;
    unsigned char *resultptr; unsigned char *resultrowptr;




    for(count = 0 , resultrowptr = result , resultptr = result + ncols * (nrows - 1); 
	count < ncols; resultptr++  , resultrowptr++  , count++ ) {
	 *resultrowptr =  *resultptr = (unsigned char)0;
    }

    for(count = 0 , resultptr = result , resultrowptr = result + ncols - 1; 
	count < nrows; count++  , resultptr += ncols , resultrowptr += ncols) {
	 *resultptr =  *resultrowptr = (unsigned char)0;
    }




    for(rowcount = 1 , magrowptr = mag + ncols + 1 , gxrowptr = gradx + ncols + 1 , 
	gyrowptr = grady + ncols + 1 , resultrowptr = result + ncols + 1; 
	rowcount <= nrows - 2; 
	rowcount++  , magrowptr += ncols , gyrowptr += ncols , gxrowptr += ncols , 
	resultrowptr += ncols) {
	for(colcount = 1 , magptr = magrowptr , gxptr = gxrowptr , gyptr = gyrowptr , 
	    resultptr = resultrowptr; colcount <= ncols - 2; 
	    colcount++  , magptr++  , gxptr++  , gyptr++  , resultptr++ ) {
	    m00 =  *magptr;
	    if (m00 == 0) {
		 *resultptr = (unsigned char)255;
	    }
	    else  {
		xperp =  -(gx =  *gxptr) / ((float)m00);
		yperp = (gy =  *gyptr) / ((float)m00);
	    }

	    if (gx >= 0) {
		if (gy >= 0) {
		    if (gx >= gy)
		    {   


			z1 =  *(magptr - 1);
			z2 =  *(magptr - ncols - 1);

			mag1 = (m00 - z1) * xperp + (z2 - z1) * yperp;


			z1 =  *(magptr + 1);
			z2 =  *(magptr + ncols + 1);

			mag2 = (m00 - z1) * xperp + (z2 - z1) * yperp;
		    }
		    else 
		    {   


			z1 =  *(magptr - ncols);
			z2 =  *(magptr - ncols - 1);

			mag1 = (z1 - z2) * xperp + (z1 - m00) * yperp;


			z1 =  *(magptr + ncols);
			z2 =  *(magptr + ncols + 1);

			mag2 = (z1 - z2) * xperp + (z1 - m00) * yperp;
		    }
		}
		else 
		{   
		    if (gx >=  -gy)
		    {   


			z1 =  *(magptr - 1);
			z2 =  *(magptr + ncols - 1);

			mag1 = (m00 - z1) * xperp + (z1 - z2) * yperp;


			z1 =  *(magptr + 1);
			z2 =  *(magptr - ncols + 1);

			mag2 = (m00 - z1) * xperp + (z1 - z2) * yperp;
		    }
		    else 
		    {   


			z1 =  *(magptr + ncols);
			z2 =  *(magptr + ncols - 1);

			mag1 = (z1 - z2) * xperp + (m00 - z1) * yperp;


			z1 =  *(magptr - ncols);
			z2 =  *(magptr - ncols + 1);

			mag2 = (z1 - z2) * xperp + (m00 - z1) * yperp;
		    }
		}
	    }
	    else 
	    {   
		if ((gy =  *gyptr) >= 0)
		{   
		    if ( -gx >= gy)
		    {   


			z1 =  *(magptr + 1);
			z2 =  *(magptr - ncols + 1);

			mag1 = (z1 - m00) * xperp + (z2 - z1) * yperp;


			z1 =  *(magptr - 1);
			z2 =  *(magptr + ncols - 1);

			mag2 = (z1 - m00) * xperp + (z2 - z1) * yperp;
		    }
		    else 
		    {   


			z1 =  *(magptr - ncols);
			z2 =  *(magptr - ncols + 1);

			mag1 = (z2 - z1) * xperp + (z1 - m00) * yperp;


			z1 =  *(magptr + ncols);
			z2 =  *(magptr + ncols - 1);

			mag2 = (z2 - z1) * xperp + (z1 - m00) * yperp;
		    }
		}
		else 
		{   
		    if ( -gx >  -gy)
		    {   


			z1 =  *(magptr + 1);
			z2 =  *(magptr + ncols + 1);

			mag1 = (z1 - m00) * xperp + (z1 - z2) * yperp;


			z1 =  *(magptr - 1);
			z2 =  *(magptr - ncols - 1);

			mag2 = (z1 - m00) * xperp + (z1 - z2) * yperp;
		    }
		    else 
		    {   


			z1 =  *(magptr + ncols);
			z2 =  *(magptr + ncols + 1);

			mag1 = (z2 - z1) * xperp + (m00 - z1) * yperp;


			z1 =  *(magptr - ncols);
			z2 =  *(magptr - ncols - 1);

			mag2 = (z2 - z1) * xperp + (m00 - z1) * yperp;
		    }
		}
	    }



	    if ((mag1 > 0.000000000000000e+00) || (mag2 > 0.000000000000000e+00))
	    {   
		 *resultptr = (unsigned char)255;
	    }
	    else 
	    {   
		if (mag2 == 0.000000000000000e+00)
		     *resultptr = (unsigned char)255;
		else 
		     *resultptr = (unsigned char)128;
	    }
	}
    }
}

#line 787 "Canny5.cc"
Platform::Platform(unsigned int _idcnt, i_img_receiver (&ImgIn), i_img_sender (&ImgOut))
    : _specc::behavior(_idcnt), ImgIn(ImgIn), ImgOut(ImgOut),
    _scc_const_port_0(2ul),
    _scc_const_port_1(2ul),
    canny(++_IDcnt, q1, q2),
    din(++_IDcnt, ImgIn, q1),
    dout(++_IDcnt, q2, ImgOut),
    q1(_scc_const_port_0),
    q2(_scc_const_port_1)
{   
}

Platform::~Platform(void)
{   
}

#line 829 "Canny5.sc"
void Platform::main()
{   
    { _specc::fork _scc_fork_0(&din), _scc_fork_1(&canny), _scc_fork_2(&dout); _specc::par(
	    &_scc_fork_0, 
	    &_scc_fork_1, 
	    &_scc_fork_2, ((_specc::fork*)0));
    }
}

#line 814 "Canny5.cc"
Main::Main(unsigned int _idcnt)
    : _specc::class_type(_idcnt),
    _scc_const_port_0(2ul),
    _scc_const_port_1(2ul),
    monitor(_IDcnt, q2),
    platform(_IDcnt, q1, q2),
    q1(_scc_const_port_0),
    q2(_scc_const_port_1),
    stimulus(_IDcnt, q1)
{   
}

Main::~Main(void)
{   
}

#line 847 "Canny5.sc"
int Main::main(void)
{   
    { _specc::fork _scc_fork_0(&stimulus), _scc_fork_1(&platform), _scc_fork_2(&monitor); _specc::par(
	    &_scc_fork_0, 
	    &_scc_fork_1, 
	    &_scc_fork_2, ((_specc::fork*)0));
    }
    return 0;
}

#line 842 "Canny5.cc"
Main _scc_main(_IDcnt);

int main(void)
{   
    int _scc_main_return;
    
    _specc::start();
    _scc_main_return = _scc_main.main();
    _specc::end();
    return(_scc_main_return);
}

void _scc_bit4_err_handle(
    const _bit4& bit4vec)
{   
    char temp_bits[1024], *p;
    p=bit2str(2,&temp_bits[1023], bit4vec);
    _specc::abort(
	"ERROR:\t Casting a bit4 vector failed \n"
	"Bit4 vector contains X/Z values %s\n"
	"Simulation aborted.\n", p);
	
}

//////////////////////////////////////////////////////////////////////
// End of file Canny5.cc
//////////////////////////////////////////////////////////////////////
