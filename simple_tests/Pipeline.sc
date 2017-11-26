//
// Pipeline.sc:
// ------------
//
// author:	Rainer Doemer
// last update:	08/02/05
//
// note:	this example simply demonstrates the execution of
//		a pipeline; it also shows that piped variables
//		are automatically buffered;


#include <sim.sh>
#include <stdio.h>
#include <stdlib.h>


behavior Stage1(out int p1, out int p2, in int p3)
{
	int	Stage1ExecCount = 0;

	void main(void)
	{
	int		t1, t2;
	sim_time_string	buf;

	Stage1ExecCount++;
	printf("Stage1 execution #%d at time %s\n",
		Stage1ExecCount, time2str(buf, now()));
	printf("Stage1 input:  p3 = %d\n", p3);
	t1 = p3 + 1000 + 10000 * Stage1ExecCount;
	t2 = p3 + 2000 + 10000 * Stage1ExecCount;
	printf("Stage1 output: p1 = %d, p2 = %d\n", t1, t2);
	waitfor(5);
	p1 = t1;
	p2 = t2;
	}
};


behavior Stage2(in int p1, out int p2)
{
	int	Stage2ExecCount = 0;

	void main(void)
	{
	int		t;
	sim_time_string	buf;

	Stage2ExecCount++;
	printf("Stage2 execution #%d at time %s\n",
		Stage2ExecCount, time2str(buf, now()));
	printf("Stage2 input:  p1 = %d\n", p1);
	t = p1 + 100000 * Stage2ExecCount;
	printf("Stage2 output: p2 = %d\n", t);
	waitfor(10);
	p2 = t;
	}
};


behavior Stage3(in int p1, in int p2, out int p3)
{
	int	Stage3ExecCount = 0;

	void main(void)
	{
	int		t;
	sim_time_string	buf;

	Stage3ExecCount++;
	printf("Stage3 execution #%d at time %s\n",
		Stage3ExecCount, time2str(buf, now()));
	printf("Stage3 input:  p1 = %d, p2 = %d\n", p1, p2);
	t = p1 + p2;
	printf("Stage3 output: p3 = %d\n", t);
	waitfor(7);
	p3 = t;
	}
};


behavior Pipeline(in int a, out int b)
{
	piped int	From1To2,
			From2To3;
	piped piped int	From1To3;

	Stage1	s1(From1To2, From1To3, a);
	Stage2	s2(From1To2, From2To3);
	Stage3	s3(From1To3, From2To3, b);

	void main(void)
	{
	pipe{	s1.main();
		s2.main();
		s3.main();
		}
	}
};


behavior TimeOut(void)
{
	void main(void)
	{
			
	waitfor(100);
	puts("TimeOut! Exiting...");
	sim_exit(0);
	}
};


behavior Main		/* Main behavior */
{
	int		i, o;
	Pipeline	p(i, o);
	TimeOut		t;

	int main(void)
	{
	puts("Pipeline:");
	i = 42;
	printf("Pipeline input:  %d\n", i);

	par {	p.main();
		t.main();
		}

	/* never reached */
	return(0);
	}
};

// EOF
