//
// Adder.sc:
// ---------
//
// author:	Jieneng Yang
// last update:	04/12/17
//


#include <stdlib.h>
#include <stdio.h>


interface IS
{
    void Send(char);
};
interface IR
{
    char Receive(void);
};


channel C implements IS, IR
{
    event Req;
    char Data;
    event Ack;
    void Send(char X)
    { 
        Data = X;
        notify Req;
        wait Ack;
    }
    char Receive(void)
    { 
        char Y;
        wait Req;
        Y = Data;
        notify Ack;
        return Y;
    }
};


behavior Prod(IS Port)
{
    char X;
    int i;
    //char S[] = "Apples and Oranges";
    char S[18] = {'A','p','p','l','e','s',' ','a','n','d',' ','O','r','a','n','g','e','s'};
    void main(void)
    {
        printf("Producer starts.\n");
        for(i = 0; i < 18; i++)
        {
        X = S[i];
        printf("Producer sends '%c'\n",X);
        Port.Send(X);
        }
        printf("Producer ends");
    }
};

behavior Cons(IR Port)
{
    char Y;
    int j;
    void main(void)
    {   
        printf("Consumer starts.\n");
        for(j = 0; j < 18; j++)
        {
        Y = Port.Receive();
        printf("Consumer received '%c'.\n",Y);
        }
        printf("Consumer starts.\n");
    }
};




behavior Main(void)
{
	C	Chan;
	Prod    S(Chan);
	Cons	R(Chan);

	int main(void)
	{
        printf("Main starts.\n");
	par {	S.main();
		R.main();
		}
	printf("Main Done.\n");
	return 0;
	}
};

// EOF

