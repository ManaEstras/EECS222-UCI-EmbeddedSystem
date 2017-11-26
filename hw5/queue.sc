// queue.sc:     simple example of using typed queue channel
// author:       Rainer Doemer
// last update:  10/09/09

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sim.sh>

#include <c_typed_queue.sh>	/* make the template available */

typedef int block[64];	/* define our communication data type */

// define the sender/receiver interfaces for the data type
DEFINE_I_TYPED_SENDER(block, block)	// creates interface i_block_sender
DEFINE_I_TYPED_RECEIVER(block, block)	// creates interface i_block_receiver
DEFINE_I_TYPED_TRANCEIVER(block, block)	// creates interface i_block_tranceiver

// define the double_handshake channel for the data type
DEFINE_C_TYPED_QUEUE(block, block)	// creates channel
					// c_block_queue

behavior Sender(i_block_sender Port)		// the sender behavior
{
    void main(void)
    {
	int Data[64] = { 42, 42, 42 };

	// ...
	Port.send(Data);
	// ...
    }
};

behavior Receiver(i_block_receiver Port)	// the receiver behavior
{
    void main(void)
    {
	int Data[64];

	// ...
	Port.receive(&Data);
	// ...
    }
};

behavior Main	// let the sender send data to the receiver through the channel
{
    c_block_queue	c(10ul);
    Receiver		r(c);
    Sender		s(c);

    int main(void)
    {
	par { r.main();
	      s.main();
	     }
	return 0;
    }
};

// EOF queue.sc
