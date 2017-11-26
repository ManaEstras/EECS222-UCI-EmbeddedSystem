//
// Timing.sc:
// ----------
//
// author:	Rainer Doemer
// last update:	08/02/05
//
// note:	this example demonstrates the use of the do-timing
//		construct for protocol modeling; the access protocol
//		to a RAM is encapsulated in an adapter channel
//		which connects a test component with the actual RAM;
//		this example also shows how to replace the default
//		range check with a user-defined one (see at the end)


#include <sim.sh>
#include <stdio.h>
#include <stdlib.h>


//#define USE_MY_RANGE_CHECKER	/* uncomment to enable our range check */


// high-level interface to the RAM

interface RAM_I
{
bit[7:0] ReadByte(
	bit[15:0]	Address);

void WriteByte(
	bit[15:0]	Address,
	bit[7:0]	Data);
};


// the RAM itself

behavior RAM(
	in bit[15:0]	Addr,
	inout bit[7:0]	Data,
	in bit[0:0]	Rd,
	in bit[0:0]	Wr,
	in event	CS)
{
bit[7:0]	Storage[1<<15];

void access_cycle(void)
	{
	wait(CS);
	if (Rd)
	   { waitfor(10);
	     Data = Storage[Addr];
	    }
	else
	   { if (Wr)
		{ waitfor(12);
		  Storage[Addr] = Data;
		 } /* fi */
	    }
	}

void main(void)
	{
	while(true)
	   { access_cycle();
	    }
	}
};


// the adapter channel

channel RAM_C(
	out bit[15:0]	ABus,
	inout bit[7:0]	DBus,
	out bit[0:0]	RMode,
	out bit[0:0]	WMode,
	out event	CSelect) implements RAM_I
{

bit[7:0] ReadByte(
	bit[15:0]	Address)
	{
	bit[7:0]	MyData;

	do	{ t1: { ABus = Address;
			}
		  t2: { RMode = 1;
			WMode = 0;
			notify CSelect;
			waitfor(11);	/* pass the range check */
//			waitfor(1);	/* fail the range check */
			}
		  t3: { }
		  t4: { MyData = DBus;
			}
		  t5: { ABus = 0;
			}
		  t6: { RMode = 0;
			WMode = 0;
			waitfor(10);
			}
		  t7: { }
		}
	timing
		{ range(t1; t2;  0;   );
		  range(t1; t3; 10; 20);
		  range(t2; t3; 10; 20);
		  range(t3; t4;  0;   );
		  range(t4; t5;  0;   );
		  range(t5; t7; 10; 20);
		  range(t6; t7;  5; 10);
		}
	return(MyData);
	}

void WriteByte(
	bit[15:0]	Address,
	bit[7:0]	Data)
	{
	do	{ t1: { ABus = Address;
			}
		  t2: { RMode = 0;
			WMode = 1;
			notify CSelect;
			waitfor(10);
			}
		  t3: { DBus = Data;
			}
		  t4: { waitfor(5);
			}
		  t5: { ABus = 0;
			}
		  t6: { RMode = 0;
			WMode = 0;
			waitfor(10);
			}
		  t7: { DBus = 0;
			}
		}
	timing
		{ range(t1; t2;  0;   );
		  range(t1; t3; 10; 20);
		  range(t2; t3; 10; 20);
		  range(t3; t4;  0;   );
		  range(t4; t5;  0;   );
		  range(t5; t7; 10; 20);
		  range(t6; t7;  5; 10);
		}
	}
};


// test vector generator and monitor

behavior Tester(RAM_I TestPort)
{
	void main(void)
	{
	int		Value;
	sim_time_string	buf;

	puts("RAM test:");

	printf("Time =%5s: writing value 42 to address 100...\n",
		time2str(buf, now()));
	TestPort.WriteByte(100, 42);

	printf("Time =%5s: writing value 27 to address 101...\n",
		time2str(buf, now()));
	TestPort.WriteByte(101, 27);

	printf("Time =%5s: reading from address 100...",
		time2str(buf, now()));
	Value = TestPort.ReadByte(100);
	printf(" value = %d\n", Value);

	printf("Time =%5s: reading from address 101...",
		time2str(buf, now()));
	Value = TestPort.ReadByte(101);
	printf(" value = %d\n", Value);

	sim_exit(0);	/* done, quit the program */
	}
};


// test bench

behavior Main
{
	bit[16]	ABus;
	bit[8]	DBus;
	bit[1]	RMode,
		WMode;
	event	CSelect;

	RAM	TestRAM(ABus, DBus, RMode, WMode, CSelect);
	RAM_C	Adapter(ABus, DBus, RMode, WMode, CSelect);
	Tester	Driver(Adapter);

	int main(void)
	{
	par {	TestRAM.main();
		Driver.main();
	     }
	return(0);
	}
};


#ifdef USE_MY_RANGE_CHECKER

// the following defines our own range checker
// to be used instead of the 'scc'-generated default
// (this one does print a message, but does not abort the program)

void _scc_range_check(
	sim_time	Time1, 
	sim_time	Time2, 
	sim_time	Min, 
	bool		CheckMin, 
	sim_time	Max, 
	bool		CheckMax, 
	const char	*Label1, 
	const char	*Label2, 
	const char	*File, 
	unsigned int	Line)
	{
	sim_time_string	buf;
	if (  ((CheckMin) && (Time2 - Time1 < Min))
	    ||((CheckMax) && (Time2 - Time1 > Max)))
	   { printf("\nERROR:\tTiming check failed in line %d\n"
			"\tin file \"%s\":\n",
			Line, File);
	     printf("\tTime stamp at label %s is %s.\n",
			Label1, time2str(buf, Time1));
	     printf("\tTime stamp at label %s is %s.\n",
			Label2, time2str(buf, Time2));
	     printf("\tMinimum difference is     %s.\n",
			time2str(buf, Min));
	     printf("\tMaximum difference is     %s.\n",
			time2str(buf, Min));
	     printf("\tActual difference is      %s.\n",
			time2str(buf, Time2 - Time1));
	     /* abort(); */
	    } /* fi */
	} /* end of _scc_range_check */

#endif /* USE_MY_RANGE_CHECKER */

// EOF
