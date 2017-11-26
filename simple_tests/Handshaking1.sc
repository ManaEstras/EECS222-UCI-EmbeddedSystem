//
// Handshaking1.sc:
// ----------------
//
// author:	Rainer Doemer
// last update:	05/18/01
//
// note:	this file is an example for a communication protocol
//		using two-way handshaking; the sender sends 10 blocks
//		of data to the receiver, performing a double-handshake
//		for each single data element transmitted;
//		this communication protocol is only safe when both
//		partners call the channel at the same time;
//		(see also Handshaking2.sc)


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

	void Send(int i)
	{
		Data = i;
		notify Valid;
		wait Ack;
	}

	int Receive(void)
	{
		int	i;

		wait Valid;
		i = Data;
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
