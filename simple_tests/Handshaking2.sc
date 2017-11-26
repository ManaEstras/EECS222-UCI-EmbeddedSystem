//
// Handshaking2.sc:
// ----------------
//
// author:	Rainer Doemer
// last update:	05/18/01
//
// note:	this file is an example for a communication protocol
//		using two-way handshaking; the sender sends 10 blocks
//		of data to the receiver, performing a double-handshake
//		for each single data element transmitted;
//		compared to Handshaking1.sc, this example uses boolean
//		flags to guard the actual handshake events; this avoids
//		potential deadlocks when using 'waitfor' since it ensures
//		that no event is lost, for example, if the receiver calls
//		the communication channel at a later time than the sender;

#include <stdio.h>

interface TwoWayHandShaking
{
	void Send(int i);
	int Receive(void);
};

channel HandShaking implements TwoWayHandShaking
{
	int	Data;
	event	Valid,
		Ack;
	bool	ValidFlag = false,
		AckFlag = false;

	void Send(int i)
	{
		Data = i;
		ValidFlag = true;
		notify Valid;
		while(! AckFlag)
			wait Ack;
		ValidFlag = false;
		notify Valid;
		while(AckFlag)
			wait Ack;
	}

	int Receive(void)
	{
		int	i;

		while(! ValidFlag)
			wait Valid;
		i = Data;
		AckFlag = true;
		notify Ack;
		while(ValidFlag)
			wait Valid;
		AckFlag = false;
		notify Ack;

		return i;
	}
};

behavior Sender(TwoWayHandShaking Port)
{
	void main(void)
	{
	int	i, j;
	int	d, s = 0;

	for(i=0; i<10; i++)
	{
		waitfor 5;
		for(j=0; j<1024; j++)
		{	d = i*100000+j;
			Port.Send(d);
			s += d;
		}
	}
	printf("Sender done. Checksum = %d\n", s);
	}
};

behavior Receiver(TwoWayHandShaking Port)
{
	void main(void)
	{
	int	i, j;
	int	d, s = 0;

	for(i=0; i<10; i++)
	{
		waitfor 10;
		for(j=0; j<1024; j++)
		{	d = Port.Receive();
			s += d;
		}
	}
	printf("Receiver done. Checksum = %d\n", s);
	}
};

behavior Main(void)
{
	HandShaking	C;
	Sender		S(C);
	Receiver	R(C);

	int main(void)
	{
	par {	S.main();
		R.main();
		}
	printf("Exiting.\n");
	return 0;
	}
};

// EOF
